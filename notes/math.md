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

### Functions and Parameters
For a function $f$ that represents a model, we write $f(x, y; a, b)$ when $x$
and $y$ are input variables to the model and $a$ and $b$ are model parameters.

### Cortical Magnification
We define the cortical magnification function $m(x, y)$, with $x$ and $y$ being
the coordinates of a point in the visual field typically in degrees of visual
angle, to be a function that yields the magnification factor, typically in
$\mbox{mm}^2$ of cortex per $\mbox{degree}^2$ of visual angle. This notably
differs from the linear cortical magnification such as that proposed by Horton
and Hoyt (1991) in their equation $m_{\sqrt{\mbox{HH}}}(r; a, b) = a / (b +
r)$ for an eccentricity $r = \sqrt{x^2 + y^2}$ and parameters $a$ in mm and $b$
in degrees. (Note that we use the identifier $m_{\mbox{HH}}$ to refer to the
areal cortical magnification model and the identifier $m_{\sqrt{\mbox{HH}}}$ to
refer to the linear version.) Such a function can be converted from a linear
cortical magnification to an (areal) cortical magnification function by
squaring it; here, we use the areal version of the Horton and Hoyt model:

$$ m_\mbox{HH}(r; a, b) = \left(\frac{a}{b+r}\right)^2 $$

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
       &=& \int_{-\pi}^{\pi} \int_{0}^{r} \rho \\, m(\rho\cos\theta, \rho\sin\theta) \\, \mathrm{d}\rho \\, \mathrm{d}\theta \\
       &=& A_0 \int_{-\pi}^{\pi} \int_{0}^{r} \rho \\, f(\rho\cos\theta, \rho\sin\theta) \\, \mathrm{d}\rho \\, \mathrm{d}\theta \\
\end{aligned} $$

Depending on the functions $m$ and $f$, this construction may not have a closed
form; though in the specific case of a radial cortical magnification function
(i.e., a cortical magnification function that depends only on the eccentricity
$r$ and not on the specific visual field position $(x,y)$ ), we can simplify
this slightly more:

$$ \begin{align}
M(r) &=& A_0 F(r) \\
     &=& \int_{-\pi}^{\pi} \int_{0}^{r} \rho \\, m(\rho) \\, \mathrm{d}\rho \\, \mathrm{d}\theta \\
     &=& 2 \pi \int_{0}^{r} \rho \\, m(\rho) \\, \mathrm{d}\rho \\
     &=& 2 \pi A_0 \int_{0}^{r} \rho \\, f(\rho) \\, \mathrm{d}\rho
\end{align} $$

### The Horton and Hoyt (1991) model
In the specific case of the model $m_\mbox{HH}(r; a, b) = (a / (b + r))^2$, we
can solve the integral directly:

$$ \begin{align}
M_\mbox{HH}(r; a, b) &=& \int_{-\pi}^{\pi} \int_{0}^{r} \rho \\, \left( \frac{a}{b+\rho} \right)^2 \\, \mathrm{d}\rho \\, \mathrm{d}\theta \\
    &=& 2 \pi a^2 \int_{0}^{r} \rho \\, (b+\rho)^{-2} \\, \mathrm{d}\rho \\
    &=& 2 \pi a^2 \left( \log\left(\frac{b+r}{b}\right) - \frac{r}{b + r} \right)
\end{align} $$

Given the closed form of this integral, we can reparameterize the model in
terms of $A_0$ by observing that $A_0 = M_\mbox{HH}(R; a, b)$ and solving for
the parameter $a$:

$$ \begin{aligned}
A_0 &=& M_{\mbox{HH}}(R; a, b) \\
  &=& 2 \pi a^2 \left( \log\left(\frac{b+R}{b}\right) - \frac{R}{b + R} \right) \\
2 \pi a^2 &=& \frac{A_0}{\log\left(\frac{b+R}{b}\right) - \frac{R}{b + R}} \\
a &=& \sqrt{ \frac{A_0}{2 \pi \left(\log\left(\frac{b+R}{b}\right) - \frac{R}{b + R}\right)} }
\end{aligned} $$

This substitution gives us a version of original model based on $A_0$ and $R$ instead of the parameter $a$:

$$ m_{\mbox{HH}}(r; b, A_0, R) = \frac{A_0}{2 \pi \left(\log\left(\frac{b+R}{b}\right) - \frac{R}{b + R}\right) (b + r)^2} $$

The Horton and Hoyt model is closely related to the [reciprocal
distribution](https://en.wikipedia.org/wiki/Reciprocal_distribution). This
distribution is characterized by the probability density function
$f_\mbox{recip}(x)$ and cumulative density function $F_\mbox{recip}(x)$, below:

$$ \begin{aligned}
f_\mbox{recip}(x; \alpha, \beta) &=& \frac{1}{x \\, \log\left(\frac{\beta}{\alpha}\right)} \\
F_\mbox{recip}(x; \alpha, \beta) &=& \frac{\log(x) - \log(\alpha)}{\log(\beta) - \log(\alpha)} 
\end{aligned} $$

The linear version of the Horton and Hoyt model is essentially a reciprocal
distribution that has been shifted along the $x$-axis:

$$ \begin{aligned}
m_{\sqrt{\mbox{HH}}}(r; a, b) &=& \sqrt{A_0} \\, f_\mbox{recip}(r + b; b, b + R) \\
   &=& \sqrt{A_0} \\, \left((r + b) \\, \log\left(\frac{b + R}{b}\right)\right)^{-1} \\
   &=& \frac{\sqrt{A_0} \log^{-1}\left(\frac{b+R}{b}\right)}{r + b}
\end{aligned} $$



## Fitting cortical magnification


