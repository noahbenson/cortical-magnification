# Notes on the Mathematics of the Cortical Magnification Function

## Notation

We define the cortical magnification function $m(x, y)$, with $x$ and $y$ being
the coordinates of a point in the visual field typically in degrees of visual
angle, to be a function that yields the magnification factor, typically in
mm$^2$ of cortex per degree$^2$ of visual angle. This notably differs from the
linear cortical magnification such as that proposed by Horton and Hoyt (1991)
in their equation $m_\mbox{HH}(r) = a / (b + r)$ for an eccentricity $r =
\sqrt{x^2 + y^2}$ and parameters $a$ in mm and $b$ in degrees. Such a function
can be converted from a linear cortical magnification to an (areal) cortical
magnification function by squaring it.

We use the variable $r$ to indicate the visual eccentricity ($r = \sqrt{x^2 +
y^2}$).

When dealing with cortical magnification functions, we assume that they are
limited to a single visual area on cortex (such as V1) that represents the
entire visual field (out to the maximum eccentricity $R$). Although this
simplifying assumption is not true (the human visual field is not a uniform
disk), the cortical magnification at high eccentricities is sufficiently small
that the effect is minimal.


## The cortical magnification function as a distribution

The cortical magnification function can be interpreted as analogous to a
probability distribution. To arrive at this construction, consider the
distrubution of pRF centers on the cortical surface. The probability that a
point randomly chosen from the cortical surface has a particular pRF center
$(x, y)$ is given by this hypothetical distribution. The probability density
function $f(x, y)$ of this distribution is simply the cortical magnification
function $m$ divided by the overall size of the represented visual area, which
we denote $A_0$:

$$ m(x, y) = A_0 f(x,y) $$.
