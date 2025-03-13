# -*- coding: utf-8 -*-
###############################################################################
# fitting.py

"""Code and functions for fitting models of Cortical Magnification.

The fitting function `fitcmag_cumarea` is the primary fitting utility, and it
requires a model function that specify's the model of cumulative area in terms
of eccentricity. (This concept of a "cumulative area" is often shortened in the
code to `cumarea`.) The cumulative area at an eccentricity `r` is the surface
area of the modeled visual area that is mapped to a visual eccentricity less
than `r` (in other words, the cortical surface area of the radius-`r` disk in
the visual field).

All models of cumulative area can be specified as continuous probability
distributions defined on the interval [0, R] where R is the maximum
eccentricity of the modeled area. (The field of view of a visual area may not
necessarily be circular, but because the cortical magnification is so low at
high eccentricities, the approximation of a disk is very close.) Probability
distributions always have an integral of 1, but a traditional cortical
magnification function (such as that of Horton & Hoyt, 1991) can be truncated
at R then normalized by the modeled visual area's total surface area to make it
compatible. For a probability distribution that meets these criteria, its
cumulative density function `F(r)` represents the normalized cumulative area:
the fraction of the visual area's total surface area that is mapped to the disk
of radius `r`. In this way, all cumulative area models handled by this module
use the form `M(r) = A * F(r)` where `A` is the total surface area of the
modeled visual area and `F(r)` is a valid cumulative distribution function on
the interval [0, R].

To transform a traditional cortical magnificaion model into a model of
cumulative area, the integral of the model over the visual field must be
found. For a function `m(r)` that, like Horton & Hoyt's model, expresses the
cortical magnification in square mm per square degree at a given eccentricity
`r` (in degrees), the cumulative area function `M(r)` is computed using the
following integral:

.. math::

    \int_{-\pi}^{\pi} \int_{0}^{r} u \, m(u) \, \mathrm{d}u \, \mathrm{d}\theta
    = 2\pi \int_{0}^{r} u \, m(u) \, \mathrm{d}u

Closed forms exist for some of the potential choices of `m` in the above
integral, but finding them is not trivial. If, instead, the cumulative area
function `M(r)` is known, then the cortical magnification function `m(r)` must
be `M'(r) / (2 pi r)`. If `M(r)` follows this module's convention (`M(r) = A *
F(r)`, see above), then the model's traditional cortical magnification function
is `m(r) = (A / (2 pi r)) * f(r)` where `f(r)` is the probability density
function of the distrubition for which `F(r)` is the cumulative density
function. Models can be encoded in this module using the `CMagRadialModel` type
in the `cmag.models` subpackage.
"""


# Fitting C.Mag. ##############################################################

def fit_cumarea(surface_areas, eccen,
                formfn, params0,
                total_area=1,
                fov=Ellipsis,
                hemifields=2,
                fit_total_area=True,
                lossfn='mse',
                weights=None,
                argtx=Ellipsis,
                **kwargs):
    """Fits a model of cortical magnification to a set of data using the method
    of cumulative surface area.

    `fitcmag_cumarea(surface_area, eccen, f, x0)` finds the parameters `x` that
    minimize the difference between `f(x)` and the cumulative distribution of
    surface area implied by `eccen` and `surface_area`. The first two
    parameters must be vectors of values corresponding to the midgray surface
    area (`surface_area`) and the eccentricity, in degrees of the visual field,
    of the pRF center of each cortical surface vertex.

    The method of cumulative area fits the cortical magnification function
    using the following steps:
     1. Sort the vertices by their eccentricity values.
     2. Calculate the cumulative sum of the surface areas of the vertices using
        the eccentricity-sorted vertex list.
     3. Find the parameters that minimize the difference between the
        eccentricity-sorted cumulative surface area across vertices and the
        predicted cumulative surface area expressed by the `formfn`.
    An implication of this method is that the `formfn` argument must be a
    function `f` such that `f(r, *params)` returns the predicted cumulative
    surface area of the portion of V1 whose eccentricity is less than `r` (in
    degrees).

    Parameters
    ----------
    surface_area : array or tensor of floats
        The midgray surface area, in square mm, of each vertex on the cortical-
        surface of the relevant visual area.
    eccen : array or tensor of floats
        The eccentricity, in degrees of the visual field, of the pRF center of
        each vertex on the cortical surface of the relevant visual area.
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
    total_area : positive float, optional
        The surface area of the predicted map. This is effectively a scale or
        gain parameter on the cumulative distribution function. The default
        value is 1, indicating that no scaling is performed. This option is
        ignored unless the `formfn` is a `cmag.models.CMagModel` object.
    fov : Ellipsis or positive float, optional
        The diameter of the field of view from which the visual area receives
        input. Note that this should be twice the maximum eccentricity from
        which the area receives input. The default value, `Ellipsis`, indicates
        that the value `cmag.hcp.config.fov` should be used. This option is
        ignored unless the `formfn` is a `cmag.models.CMagModel` object.
    hemifields : float, optional
        The number of hemifields included in the calculation of `surface_area`,
        with the default being 2. A value of 1 indicates that only one
        hemifield / hemisphere was be included, thus `total_area` is actually
        the surface area of only a LH or only a RH visual area. Wedges can also
        be specified (e.g., `hemifields=0.5` is equivalent to one quadrant of
        the visual field). This option is ignored unless the `formfn` is a
        `cmag.models.CMagModel` object.
    fit_total_area : boolean, optional
        Whether to fit the total area of the visal area (i.e., the cumulative
        radial area at an eccentricity of `fov / 2`) as a parameter. If this
        option is set to `True` (the default), then `total_area` is used as an
        initial estimate, and the actual total area is fit as part of the 
        optimization; an attribute `total_area` is added to the returned report
        object with the fitted value. If this option is set to `False`, then
        the `total_area` is taken as absolute and not fitted.
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
    from numpy import asarray, argsort, sum, cumsum, mean, sqrt, square
    from scipy.optimize import minimize
    from .models import CMagModel, CMagRadialModel
    sarea = asarray(surface_areas)
    eccen = asarray(eccen)
    params0 = asarray(params0, dtype=eccen.dtype)
    ii = argsort(eccen)
    sarea = sarea[ii]
    eccen = eccen[ii]
    cumsa = cumsum(sarea)
    if argtx is None:
        argtx = (lambda a:a, lambda a:a)
    elif argtx is Ellipsis:
        if hasattr(formfn, 'argtx'):
            argtx = formfn.argtx
    else:
        params0 = argtx[0](params0)
    if lossfn == 'rss':
        def lossfn(gold, pred):
            return sum((gold - pred)**2)
    elif lossfn == 'mse':
        if weights is None:
            def lossfn(gold, pred):
                return mean((gold - pred)**2)
        else:
            wsum = sum(weights)
            def lossfn(gold, pred):
                return sum(weights * (gold - pred)**2) / wsum
    if issubclass(formfn, CMagModel):
        import torch
        ecctns = torch.as_tensor(eccen)
        if issubclass(formfn, CMagRadialModel):
            fn = formfn.radial_cumarea
        else:
            from warnings import warn
            warn(
                f"cortical magnification model of type {type(formfn)} is a"
                f" CMagModel but not a CMagRadialModel; using model(x=eccen,"
                f" y=0)")
            y = ecctns * 0
            fn = lambda x, *args: formfn.areal_cumarea(x, y, *args)
        def stepfn(params, *args):
            if fit_total_area:
                txparams = torch.as_tensor(argtx[1](params[:-1]))
                totarea = torch.square(torch.as_tensor(params[-1]))
            else:
                txparams = torch.as_tensor(argtx[1](params))
                totarea = total_area
            pred = fn(ecctns, totarea, fov, hemifields, *txparams, *args)
            l = lossfn(cumsa, pred.numpy())
            return l
        if fit_total_area:
            params0 = list(params0)
            params0.append(sqrt(total_area))
    elif fit_total_area:
        raise ValueError(
            "total_area='fit' can only be used with CMagModel forms")
    else:
        def stepfn(params, *args):
            txparams = argtx[1](params)
            pred = formfn(eccen, *txparams, *args)
            return lossfn(cumsa, pred)
    r = minimize(stepfn, params0, **kwargs)
    if fit_total_area:
        total_area = square(r.x[-1])
        r.x = r.x[:-1]
        r.total_area = total_area
    r.x = argtx[1](r.x)
    return r
