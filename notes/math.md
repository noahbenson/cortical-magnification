# Notes on the Mathematics of the Cortical Magnification Function

## Notation

### Angles and Degrees
Degrees of the visual field and degrees of rotation (or degrees of polar angle)
can be ambiguous. To avoid this ambiguouty in this document, we always use
degrees of visual angle when discussing distances between points in the visual
field or eccentricity, and we always use radians when discussing polar angle
(which is often represented as degrees of rotation around the fovea). In other
words, degrees always refers to degrees of visual angle and radians always
refers to rotation around the visual field. A polar angle of 0 is considered to
be the right horizontal meridian, and a polar angle of $\pi / 2$ is considered
to be the upper vertical meridian.

### Cortical Magnification
We define the cortical magnification function $m(x, y)$, with $x$ and $y$ being
the coordinates of a point in the visual field typically in degrees of visual
angle, to be a function that yields the magnification factor, typically in
$\mbox{mm}^2$ of cortex per $\mbox{degree}^2$ of visual angle. This notably
differs from the linear cortical magnification such as that proposed by Horton
and Hoyt (1991) in their equation $m(r) = a / (b + r)$ for an eccentricity $r =
\sqrt{x^2 + y^2}$ and parameters $a$ in mm and $b$ in degrees. Such a function
can be converted from a linear cortical magnification to an (areal) cortical
magnification function by squaring it; here, we use the areal version of the
Horton and Hoyt model:

$$ m_\mbox{HH}(r) = \left(\frac{a}{b+r}\right)^2 $$

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

$$ m(x, y) = A_0 f(x,y) $$

One small issue with the above construction is that the cumulative density
function is difficult to derive as it requires integration:

$$ \begin{aligned}
M(x,y) &=& A_0 F(x,y) \\
       &=& \int \int m(x,y) \\, \mathrm{d}x \\, \mathrm{d}y \\
       &=& \int_{-\pi}^{\pi} \int_{0}^{R} r \\, m(r\cos\theta, r\sin\theta) \\, \mathrm{d}r \\, \mathrm{d}\theta \\
       &=& A_0 \int_{-\pi}^{\pi} \int_{0}^{R} r \\, f(r\cos\theta, r\sin\theta) \\, \mathrm{d}r \\, \mathrm{d}\theta \\
\end{aligned} $$

Depending on the functions $m$ and $f$, this construction may not have a closed
form; though in the specific case of a radial cortical magnification function
(i.e., a cortical magnification function that depends only on the eccentricity
$r$ and not on the specific visual field position $(x,y)$), we can simplify
this slightly more:

$$ \begin{align}
M(r) &=& A_0 F(r) \\
     &=& \int_{-\pi}^{\pi} \int_{0}^{R} r \\, m(r) \\, \mathrm{d}r \\, \mathrm{d}\theta \\
     &=& 2 \pi A_0 r \int_{0}^{R} m(r) \\, \mathrm{d}r
\end{align} $$


## Fitting cortical magnification
