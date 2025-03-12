# -*- coding: utf-8 -*-
###############################################################################
# hh91.py

"""Functions related to Horton & Hoyt's (1991) model of cortical magnification.

This namespace contains functions to calculate Horton & Hoyt's model as well as
its integral and other related measures.

For more details, see the original publication: 

Horton JC, Hoyt WF (1991) The representation of the visual field in human
  striate cortex. A revision of the classic Holmes map. Arch Ophthalmol.
  109(6):816-24. doi:10.1001/archopht.1991.01080060080030.
"""

def HH91(r, a=17.3, b=0.75, output='areal'):
    """Computes the V1 cortical magnification at the given eccentricity.

    The Horton and Hoyt (1991) model of cortical magnification predicts that
    the magnification (in square mm per square degree) at an eccentricity `r`
    is equal to `(a / (b + r))**2` where `a` was measured to be 17.3 mm and `b`
    was measured to be 0.75 degrees. (The square root of this quantity, `a / (b
    + r)`, is often called the linear cortical magnification while this
    quantity itself is sometimes called the areal cortical magnification.)

    Parameters
    ----------
    r : number or array or tensor
        The eccentricity, in degrees, at which to compute the cortical
        magnification.
    a : float, optional
        The parameter `a` of the model, in mm; the default value is 17.3, which
        was measured by Horton & Hoyt (1991).
    b : float, optional
        The parameter `b` of the model, in degrees of the visual field; the
        default value is 0.75, which was measured by Horton & Hoyt (1991).
    output : 'areal' or 'linear'
        Whether to return the areal or linear cortical magnification. The
        default is to return the areal magnification, which is just the square
        of the linear magnification.

    Returns
    -------
    float or array or tensor
        The predicted V1 cortical magnification at the given eccentricity
        value(s). If the `output` parameter is set to `'areal'`, then the units
        of the return value are square-mm / square-deg; otherwise, they are
        mm / deg.

    References
    ----------
    Horton JC, Hoyt WF (1991) The representation of the visual field in human
      striate cortex. A revision of the classic Holmes map. Arch Ophthalmol.
      109(6):816-24. doi:10.1001/archopht.1991.01080060080030.
    """
    lin_cmag = (a / (r + b))
    if output == 'areal':
        return lin_cmag ** 2
    elif output == 'linear':
        return lin_cmag
    else:
        raise ValueError(
            f"invalid output parameter; valid choices are 'areal' or 'linear'")

def HH91_integral(ecc, maxecc=None, /, a=17.3, b=0.75, hemifields=2):
    """Computes the integral of the V1 cortical magnification between two
    eccentricity values.

    The Horton and Hoyt (1991) model of cortical magnification predicts that
    the magnification (in square mm per square degree) at an eccentricity `r`
    is equal to `(a / (b + r))**2` where `a` was measured to be 17.3 mm and `b`
    was measured to be 0.75 degrees. This quantity describes the amount of the
    V1 cortical surface that is devoted to processing a square degree of the
    visual field at the given eccentricity. If we integrate this over the
    entire visual field, we get a prediction of the total surface area of V1.

    This function returns the integral from some minimum eccentricity (usually
    0) to some maximum eccentricity. This integral represents the surface area
    of V1 devoted to processing the visual field inputs between the given
    eccentricities.

    `HH91_integral(r)` returns the predicted V1 cortical surface area devoted
    to processing the central `r` degrees of the visual field.

    `HH91_integral(a, b)` returns the predicted V1 cortical surface area
    devoted to processing the ring of the visual field between `a` degrees and
    `b` degrees of eccentricity.

    See the accompanying notebook for the integration.

    Parameters
    ----------
    ecc : number or array or tensor
        The eccentricity, in degrees, at which to compute the cortical
        magnification integral. If `ecc` is given and the second argument
        `maxecc` is not given or is `None`, then `ecc` is used as the maximum
        eccentricity and 0 is used as the minimum eccentricity. Otherwise,
        `ecc` is used as the minimum eccentricity and `maxecc` is used as the
        maximum.
    maxecc : number or array or tensor, optional
        If given, then the maximum eccentricity up to which the integral is
        computed.  magnification integral. If `maxecc` is not given or is
        `None` then the first argument `ecc` is used as the maximum
        eccentricity and 0 is used for the minimum eccentricity.
    a : float, optional
        The parameter `a` of the model, in mm; the default value is 17.3, which
        was measured by Horton & Hoyt (1991).
    b : float, optional
        The parameter `b` of the model, in degrees of the visual field; the
        default value is 0.75, which was measured by Horton & Hoyt (1991).
    hemifields : float, optional
        The number of hemifields over which to integrate, with the default
        being 2. A value of 1 indicates that only one hemifield / hemisphere is
        to be included. Wedges can also be specified (e.g., `hemifields=0.5` is
        equivalent to one quadrant of the visual field).

    Returns
    -------
    float or array or tensor
        The predicted V1 cortical surface area between the given eccentricity
        values. The units are square-mm.

    References
    ----------
    Horton JC, Hoyt WF (1991) The representation of the visual field in human
      striate cortex. A revision of the classic Holmes map. Arch Ophthalmol.
      109(6):816-24. doi:10.1001/archopht.1991.01080060080030.
    """
    from torch import is_tensor
    from numpy import pi
    # Process the arguments:
    (r0,r1) = (0, ecc) if maxecc is None else (ecc, maxecc)
    # Figure out which backend we're using:
    if any(map(is_tensor, (r0, r1, b))):
        from torch import log
    else:
        from numpy import log
    # The Integral; according to Mathematica, the following input:
    #   Assuming[
    #     And[
    #       Element[{r, a, b, r0, r1}, Reals],
    #       r >= 0, r0 >= 0, r1 >= r0,
    #       a >= 0, b > 0],
    #     FullSimplify@Integrate[
    #       r*HH91[r, a, b],
    #       {th, -Pi, Pi},
    #       {r, r0, r1}]]
    # yields:
    #   2 a^2 Pi (b (-(1/(b + r0)) + 1/(b + r1)) - Log[b + r0] + Log[b + r1]).
    # Or, in more python-ic notation and simplified:
    #   2*pi * (a**2) * (b * (1/(b+r1) - 1/(b+r0)) + log((b+r1)/(b+r0)).
    # In the case of a minimum eccentricity of 0, this simplifies to:
    #   2*pi * (a**2) * (b * (1/(b+r) - 1/(b+0)) + log((b+r)/(b+0))
    #   = 2*pi * (a**2) * (b * (b - (b+r))/(b*(b+r)) + log((b+r)/b))
    #   = 2*pi * (a**2) * ((-r)/(b+r) + log((b+r)/b))
    #   = 2*pi * (a**2) * (log((b+r)/b) - r/(b+r)).
    if r0 == 0:
        b_r1 = r1 + b
        return hemifields * pi * a**2 * (log(b_r1 / b) - r1 / b_r1)
    else:
        b_r0 = b + r0
        b_r1 = b + r1
        return hemifields * pi * a**2 * (b/b_r1 - b/b_r0 + log(b_r1/b_r0))

def HH91_find_a(surfarea, ecc, max_eccen=None, /, b=0.75, hemifields=2):
    """Calculate the parameter `a` for the Horton & Hoyt (1991) model based on
    the surface area of V1 up to a certain eccentricity.

    `HH91_find_a(A, M)` finds the `a` value that corresponds to a V1 whose most
    foveal `M` degrees of eccentricity have a surface area of `A` (in square
    mm).

    `HH91_find_a(A, m, M)` finds the `a` value that corresponds to a V1 whose
    surface area between minimum eccentricity `m` and maximum eccentricity `M`
    is `A` (in square mm).

    Parameters
    ----------
    surfarea : float or array or tensor
        The surface area of the portion of V1 limited by the eccentricity
        values of the remaining arguments.
    ecc : number or array or tensor
        The eccentricity, in degrees, that limits the surface area. If `ecc` is
        given and the second argument `maxecc` is not given or is `None`, then
        `ecc` is used as the maximum eccentricity and 0 is used as the minimum
        eccentricity. Otherwise, `ecc` is used as the minimum eccentricity and
        `maxecc` is used as the maximum.
    maxecc : number or array or tensor, optional
        If given, then the maximum eccentricity limiting the surface area. If
        `maxecc` is not given or is `None` then the first argument `ecc` is
        used as the maximum eccentricity and 0 is used for the minimum
        eccentricity.
    b : float, optional
        The parameter `b` of the model, in degrees of the visual field; the
        default value is 0.75, which was measured by Horton & Hoyt (1991).
    hemifields : float, optional
        The number of hemifields over which to integrate, with the default
        being 2. A value of 1 indicates that only one hemifield / hemisphere is
        to be included. Wedges can also be specified (e.g., `hemifields=0.5` is
        equivalent to one quadrant of the visual field).

    Returns
    -------
    float or array or tensor
        The value of `a` that, combined with the parameter `b`, predicts that
        the surface area of V1 between the minimum and maximum eccentricities
        is equal to the given surface area parameter `surfarea`.

    References
    ----------
    Horton JC, Hoyt WF (1991) The representation of the visual field in human
      striate cortex. A revision of the classic Holmes map. Arch Ophthalmol.
      109(6):816-24. doi:10.1001/archopht.1991.01080060080030.

    """
    from torch import is_tensor
    from numpy import pi
    # Process the arguments:
    (r0,r1) = (0, ecc) if max_eccen is None else (ecc, max_eccen)
    # Figure out which backend we're using:
    if any(map(is_tensor, (surfarea, r0, r1, b))):
        from torch import (log, sqrt)
    else:
        from numpy import (log, sqrt)
    # The HH91_integral function (above) gives the formula for the surface
    # area, A:
    #   A = 2*pi * (a**2) * (b * (1/(b+r1) - 1/(b+r0)) + log((b+r1)/(b+r0))).
    # We can solve this for the parameter a:
    #  A / (2*pi * (b * (1/(b+r1) - 1/(b+r0)) + log((b+r1)/(b+r0)))) = a**2
    #  a = sqrt(A / (2*pi * (b * (1/(b+r1) - 1/(b+r0)) + log((b+r1)/(b+r0))))).
    b_r0 = b + r0
    b_r1 = b + r1
    denom = (hemifields * pi * (b*(1/b_r1 - 1/b_r0) + log(b_r1/b_r0)))
    return sqrt(surfarea / denom)
