# -*- coding: utf-8 -*-
###############################################################################
# models.py

"""Models of cortical magnification.
"""

import torch
import numpy as np


# CMagModel ###################################################################

class CMagModel:
    """A model of cortical magnification.
    """
    @classmethod
    def areal_cmag(cls, x, y, total_area, fov, *args):
        raise NotImplementedError
    @classmethod
    def linear_cmag(cls, x, y, total_area, fov, *args):
        m = cls.areal_cmag(x, y, total_area, fov, *args)
        return torch.sqrt(m)
    def __new__(cls, x, y, /, *args,
                total_area=1, fov=Ellipsis, hemifields=2, device=None,
                form='areal'):
        torch_inputs = torch.is_tensor(x) or torch.is_tensor(y)
        total_area = float(total_area)
        x = torch.as_tensor(x, device=device)
        y = torch.as_tensor(y, device=device)
        tmp = []
        for arg in args:
            try:
                tmp.append(torch.as_tensor(arg, device=device))
            except Exception:
                tmp.append(arg)
        if form == 'areal':
            res = cls.areal_cmag(x, y, total_area, fov, hemifields, *args)
        elif form == 'linear':
            res = cls.linear_cmag(x, y, total_area, fov, hemifields, *args)
        else:
            raise ValueError(
                "CMagModel form argument must be 'areal' or 'linear'")
        if not torch_inputs:
            res = res.detach().numpy()
        return res
class CMagRadialModel(CMagModel):
    @classmethod
    def radial_cumarea(cls, r, total_area, fov, hemifields, *args):
        raise NotImplementedError
    @classmethod
    def radial_cmag(cls, r, total_area, fov, hemifields, *args):
        raise NotImplementedError
    @classmethod
    def radial_area(cls, r, total_area, fov, hemifields, *args):
        cmag = cls.radial_cmag(
            r, total_area, fov, *args,
            hemifields=2, **kw)
        return hemifields * torch.pi * r * cmag
    @classmethod
    def areal_cmag(cls, x, y, total_area, fov, hemifields, *args):
        r = torch.hypot(x, y)
        return cls.radial_cmag(r, total_area, fov, hemifields, *args)

class hh91(CMagRadialModel):
    """A cortical magnification model based on Horton and Hoyt's (1991) model.

    The `HH91_form` function is based on Horton and Hoyt's model model `m(r) =
    (a / (b + r))**2` where `r` is the eccentricity in degrees, `a` is a scale
    parameter related to the visual area size in mm, and `b` is a shape
    parameter. The `HH91_form` function uses slightly different parameters:
    `shape`, which corresponds to the parameter `b`, `total_area`, which
    corresponds to the total surface area of the visual area, and `fov`, which
    corresponds to the diameter of the field of view from which the visual area
    receives in put (twice the maximum eccentricity).

    See also: `cmag.hh91.HH91`, `cmag.hh91.HH91_integral`,
    `cmag.hh91.HH91_find_a`.

    Parameters
    ----------
    r : float or array
        The eccentricity or eccentricities at which to calculate the cortical
        magnification.
    shape : positive float
        The shape parameter of the Horton and Hoyt (1991) model.
    total_area : positive float, optional
        The total surface area of the predicted map. The default value is 1,
        indicating that the output will be normalized as if it were the
        cumulative distribution function of a propbability distribution.
    fov : Ellipsis or positive float, optional
        The diameter of the field of view from which the visual area receives
        input. Note that this should be twice the maximum eccentricity from
        which the area receives input. The default value, `Ellipsis`, indicates
        that the value `cmag.hcp.config.fov` should be used.
    bilateral : boolean, optional
        Whether the given surface area option `total_area` represents the
        bilateral surface area or not. If `bilateral` is set to `False`, then
        it must be the surface area for a single hemisphere. The default is
        `True`.
    """
    @classmethod
    def radial_cumarea(cls, r, total_area, fov, hemifields, b=0.75):
        from .hh91 import HH91_integral, HH91_find_a
        if fov is Ellipsis:
            from .hcp.config import fov
        max_eccen = float(fov) / 2
        a = HH91_find_a(total_area, 0, max_eccen, b=b, hemifields=hemifields)
        return HH91_integral(0, r, a=a, b=b, hemifields=hemifields)
    @classmethod
    def radial_cmag(cls, r, total_area, fov, hemifields, b=0.75):
        from .hh91 import HH91, HH91_find_a
        if fov is Ellipsis:
            from .hcp.config import fov
        max_eccen = float(fov) / 2
        a = HH91_find_a(total_area, 0, max_eccen, b=b, hemifields=hemifields)
        return HH91(r, a, b)
    argtx = (np.log, np.exp)

class beta(CMagRadialModel):
    """A cortical magnification model based on the beta distribution.

    The `beta_form` of cortical magnification uses four parameters: `a` and
    `b`, which correspond to the alpha and beta parameters of the beta
    distribution; `total_area`, which corresponds to the total predicted
    surface area of the visual area; and `fov`, which corresponds to the
    maximum field of view from which the visual area receives input. Note that
    the `fov` parameter should be twice the maximum eccentricity for the area.

    Parameters
    ----------
    r : float or array
        The eccentricity or eccentricities at which to calculate the cortical
        magnification.
    a : positive float
        The shape parameter alpha of the beta distribution used to model the
        cortical magnification.
    b : positive float
        The shape parameter beta of the beta distribution used to model the
        cortical magnification.
    total_area : positive float, optional
        The surface area of the predicted map. This is effectively a scale or
        gain parameter on the cumulative distribution function. The default
        value is 1, indicating that no scaling is performed.
    fov : Ellipsis or positive float, optional
        The diameter of the field of view from which the visual area receives
        input. Note that this should be twice the maximum eccentricity from
        which the area receives input. The default value, `Ellipsis`, indicates
        that the value `cmag.hcp.config.fov` should be used.
    bilateral : boolean, optional
        Whether the given surface area option `total_area` represents the
        bilateral surface area or not. If `bilateral` is set to `False`, then
        it must be the surface area for a single hemisphere. The default is
        `True`.

    Returns
    -------
    float or array
        The predicted cortical magnification at the eccentricity values given
        in `r` according to the provided parameters.
    """
    @classmethod
    def radial_cumarea(cls, r, total_area, fov, hemifields,
                       a=2.0, b=3.0):
        from scipy.stats import beta
        if fov is Ellipsis:
            from .hcp.config import fov
        max_eccen = float(fov) / 2
        b = total_area * beta.cdf((r / max_eccen).numpy(), a, b)
        return torch.as_tensor(b)
    @classmethod
    def radial_area(cls, r, total_area, fov, hemifields, a=2.0, b=3.0):
        from torch.distributions.beta import Beta
        if fov is Ellipsis:
            from .hcp.config import fov
        max_eccen = float(fov) / 2
        return total_area * Beta(a, b).log_prob(r / max_eccen).exp()
    @staticmethod
    def _beta(a, b):
        la = torch.lgamma(a)
        lb = torch.lgamma(b)
        lab = torch.lgamma(a + b)
        return torch.exp(la + lb - lab)
    @classmethod
    def radial_cmag(cls, r, total_area, fov, hemifields, a=1.0, b=3.0):
        if fov is Ellipsis:
            from .hcp.config import fov
        max_eccen = float(fov) / 2
        const = total_area / (hemifields * torch.pi * max_eccen)
        pdf = cls.radial_area(r, 1, fov, hemifields, a=a, b=b)
        pdf_adj = pdf / r
        if a < 2:
            pdf_adj[r == 0] = torch.inf
        elif a > 2:
            pdf_adj[r == 0] = 0
        else:
            pdf_adj[r == 0] = torch.pow(max_eccen, 1 - a) / beta._beta(a,b)
        return const * pdf_adj
    argtx = (np.log, np.exp)
