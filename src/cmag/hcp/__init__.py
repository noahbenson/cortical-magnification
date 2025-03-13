# -*- coding: utf-8 -*-
###############################################################################
# hcp/__init__.py

"""Analysis code for examining the cortical magnification of the Human
Connectome Project.

The subpackage `cmag.hcp.config` can be edited to configure the paths for
loading and caching data. Most loading of HCP data is handled by the
neuropythy library, which must be configured so that it knows where to find
those data; see https://github.com/noahbenson/neuropythy.

The subpackage `cmag.hcp.data` contains functions for loading the simplified and
cached HCP data.
"""


from . import config
from . import data

# We can extract a subject list using neuropythy:
import neuropythy as ny, numpy as np
sids = np.array(ny.data['hcp_lines'].subject_list)
sids.setflags(write=False)
# To keep the data namespace clean, we delete the library imports.
del ny, np

def fit_cmag_data(data, formfn, params0,
                  fov=Ellipsis,
                  hemifields=1,
                  lossfn='mse',
                  weights=None,
                  argtx=Ellipsis,
                  filter=None,
                  labels=None):
    """Fits the cortical magnification function to an HCP hemisphere's data.

    This function is a wrapper around the `cmag.fitting.fit_cumarea` function
    that simplifies calling it on the data from a hemisphere of an HCP subject,
    as loaded using the `cmag.hcp.data.load` function. The first argument,
    `data`, should be a dictionary as loaded by this function. The remaining
    keyword arguments are identical to those of `fit_cumarea` with the
    exception of `filter` and `labels`, which lets one choose only certain
    vertices to use in fitting.

    Parameters
    ----------
    data : dict
        A dictionary of data for a particular hemisphere, as loaded by the
        `cmag.hcp.data.load` function.
    formfn : function or CMagModel
        The form of the cumulative model of cortical magnification. The
        function call `formfn(r, *params)` must return the predicted surface
        area of the portion of the relevant visual area with an eccentricity
        less than `r`. This is the cumulative version of the cortical
        magnification, not the traditional version describing the magnification
        itself. Alternatively, a `cmag.models.CMagModel` can be given, which
        allows the various optional arguments (`total_area`, `fov`, and
        `hemifields`) to be interpreted correctly (these are ignored if
        `formfn` is not a `CMagModel` object).
    params0 : sequence
        A sequence of parameter values for `formfn` that constitute the
        starting position of the minimization search. The parameters are passed
        to `formfn` along with a vector of eccentricity values `r` using the
        syntax `formfn(r, *params)`.
    hemifields : float, optional
        The number of hemifields included in the calculation of `surface_area`,
        with the default being 1. A value of 1 indicates that only one
        hemifield / hemisphere was be included, thus `total_area` is actually
        the surface area of only a LH or only a RH visual area. Wedges can also
        be specified (e.g., `hemifields=0.5` is equivalent to one quadrant of
        the visual field). This option is ignored unless the `formfn` is a
        `cmag.models.CMagModel` object.
    lossfn : function, optional
        The loss function used in the minimization. The function must accept
        two arguments: `lossfn(meas, pred)` where `meas` is the measured
        surface area provided by the parameters `prf_x`, `prf_y`, and
        `surface_area`. If not given or if the value is `None`, then the mean
        squared error (MSE) is used.
    weights : None or vector of floats, optional
        The weights to apply to the calculation of the loss function. If `None`
        or not provided, the no weights are used (i.e., all vertices are
        weighted equally). Otherwise, the weights are used by the loss function
        to calculate the weighted MSE instead of the standard MSE. If the an
        explicit function is given for `lossfn`, then this argument is ignored.
    argtx : 2-tuple of functions, optional
        A 2-tuple of `(argin, argex)` functions. The `in_params =
        argin(params)` function call must return a transformed version of
        `params` that is appropriate for use in minimization while `params =
        argex(in_params)` must return the true value of the transformed
        parameters. This argument is typically used to prevent parameter values
        from becoming negative. For example, when fitting the traditional
        Horton and Hoyt (1991) model of cortical magnification, the functional
        form is `f(r, a, b) = (a / (b + r))**2`. If `b` is non-positive, this
        can be undefined at `r = 0`, so the `argtx=(np.log, np.exp)` argument
        is used. Because `exp(in_r)` is always positive for real-valued `r`,
        this prevents the parameter used in the `formfn` from becoming invalid
        during the search.
    filter : function or None, optional
        A function that should be applied to the given data in order to obtain
        a boolean mask of values to include in the fit. For example, to include
        only vertices whose eccentricity is less than 10, one could use
        `filter=lambda data: data['eccentricity'] < 10`.
    labels : int or list of ints or None, optional
        The labels to include in the return value. If an integer label is given
        then the return value is the `scipy.optimize.OptimizeResult` object for
        that label. If a sequence of integer labels is given, then a tuple of
        `OptimizeResult` objects is returned, one per label. The default value
        is `None`, indicating that `arange(1,11)` should be used.
    **kwargs
        Any additional options are passed to the `scipy.optimize.minimize`
        function directly.

    Returns
    -------
    res : OptimizeResult
        The result of the optimization, as returned by the function
        `scipy.optimize.minimize`. The attribute `total_area` is added to the
        returned object with either the given or fitted total area, in square
        mm.

    """
    import numpy as np
    from ..fitting import fit_cumarea
    kw = dict(
        lossfn=lossfn,
        argtx=argtx,
        hemifields=hemifields)
    if not filter:
        filter = lambda dat: True
    if labels is None:
        labels = range(1,11)
    # If we were passed a 2-tuple of (lh,rh) data, join them.
    if isinstance(data, tuple) and len(data) == 2:
        from .data import joinhemis
        data = joinhemis(data)
    res = []
    for lbl in ([labels] if isinstance(labels, int) else labels):
        ii = filter(data) & (data['label'] == lbl)
        if np.sum(ii) < 5:
            r = None
        else:
            weights = None if not weights else data[weights][ii]
            sarea = data['surface_area'][ii]
            eccen = data['eccentricity'][ii]
            r = fit_cumarea(
                sarea,
                eccen,
                formfn, params0,
                total_area=np.sum(sarea) * 2,
                fov=fov,
                fit_total_area=True,
                weights=weights,
                **kw)
        res.append(r)
    return res[0] if isinstance(labels, int) else tuple(res)

__all__ = (
    "config",
    "data",
    "sids")
