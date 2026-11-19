#set page(paper: "a4", margin: 2.2cm)
#set text(size: 11pt)
#set heading(numbering: "1.")
#set math.equation(numbering: none)

#align(center)[
  #text(size: 18pt, weight: "bold")[Vibration of a Bending Beam]\
  Yizhou Lu
]

= Introduction

A beam is clamped at one end, and a force $-F$ is applied at the free end. The beam is released and vibrates.

Initial equilibrium assumptions:
- Ignore gravity.
- Small angle / deflection.

= Governing Equations

== Bending Moment

Neutral axis length: $L_"neutral axis" = R theta$

Strain is given by
$ epsilon(y) = ((y + R) theta - R theta)/(R theta) = y / R $

Stress is given by
$ sigma(y) = E epsilon(y) = E y / R $

Bending moment is the torque of the stress:
$ upright(d) M = y upright(d) F = y sigma(y) upright(d) A = E y^2 / R upright(d) A quad => quad M = integral_A E y^2 / R upright(d) A $

Area moment of inertia:
$ I = integral_A y^2 upright(d) A $

Curvature $kappa = 1/R$, so
$ M = E I kappa $

== Internal Node

Ignoring rotation, equilibrium of torque:
$ V(x) upright(d) x + M(x + upright(d) x) = M(x) $

Newton's second law (mass of the segment is $rho A upright(d) x$, acceleration is $partial^2 H / (partial t^2)$):
$ rho A upright(d) x (partial^2 H)/(partial t^2) = V(x) - V(x + upright(d) x) $

This gives
$ cases(
  V(x) = -(partial M)/(partial x),
  rho A (partial^2 H)/(partial t^2) = -(partial V)/(partial x)
) quad => quad rho A (partial^2 H)/(partial t^2) = (partial^2 M)/(partial t^2) $

$ -(partial^2)/(partial x^2)(E I kappa) = rho A (partial^2 H)/(partial t^2) $

$ E I (partial^2)/(partial x^2) (-(partial^2 H)/(partial x^2)) = rho A (partial^2 H)/(partial t^2) $

Hence we finally get the *Euler-Bernoulli Beam Equation*:

#align(center)[
#box(stroke: 0.8pt, inset: 8pt)[
  $ E I (partial^4 H)/(partial x^4) + rho A (partial^2 H)/(partial t^2) = 0 $
]
]

Curvature (small angle):
$ kappa = 1/R = (upright(d) theta)/(upright(d) s) approx (upright(d))/(upright(d) x)(-(upright(d) H)/(upright(d) x)) = -(upright(d)^2 H)/(upright(d) x^2) $

= Numerical Methods

$ cases(
  (partial U)/(partial t) = -(E I)/(rho A) (partial^4 H)/(partial x^4),
  (partial H)/(partial t) = U quad "(velocity)"
) $

Finite differences:
$ (upright(d)^2 y)/(upright(d) x^2) = (y(x + Delta x) - 2 y(x) + y(x - Delta x))/(Delta x^2) $

$ (upright(d)^3 y)/(upright(d) x^3) = 1/(Delta x^3) [ 1/2 y(x + 2 Delta x) - y(x + Delta x) + y(x - Delta x) - 1/2 y(x - 2 Delta x) ] $

$ (upright(d)^4 y)/(upright(d) x^4) = 1/(Delta x^4) [ y(x + 2 Delta x) - 4 y(x + Delta x) + 6 y(x) - 4 y(x - Delta x) + y(x - 2 Delta x) ] $

== Boundary Conditions ($t >= 0$)

Clamped end:
- Fixed position: $H|_(x = 0) = 0$ #h(1em) (1)
- No rotation: $(partial H)/(partial x)|_(x = 0) = 0$ #h(1em) (2)

Free end:
- $M(L) prop (partial^2 H)/(partial x^2)|_(x = L) = 0$ #h(1em) (3)
- $V(L) prop (partial^3 H)/(partial x^3)|_(x = L) = 0$ #h(1em) (4)

== Initial Conditions ($t = 0$)

Torque balance keeps unchanged:
$ V(x) upright(d) x + M(x + upright(d) x) = M(x) $

Equilibrium of forces:
$ V(x) = V(x + upright(d) x) = F $

Therefore
$ F + (partial M)/(partial x) = 0 $

Recall $M = -E I (partial^2 H)/(partial x^2)$. Then
$ F = E I (partial^3 H)/(partial x^3) quad => quad (partial^3 H)/(partial x^3) = F / (E I) $

Integrating:
$ (partial^3 H)/(partial x^3) = F / (E I) x + C_1 $

$ (partial H)/(partial x) = F / (2 E I) x^2 + C_1 x + C_2 $

$ H = F / (6 E I) x^3 + C_1 / 2 x^2 + C_2 x + C_3 $

Apply boundary conditions:
$ C_1 = -F / (E I) L, quad C_2 = 0, quad C_3 = 0 $

$ H = F / (6 E I) (-x^3 + 3 L x^2) $

== Ghost Nodes

Free end, condition (3):
$ (partial^2 H)/(partial x^2)|_(x = L) = 0 => H(L + Delta x) - 2 H(L) + H(L - Delta x) = 0 $
$ => H(L + Delta x) = 2 H(L) - H(L - Delta x) $

Free end, condition (4):
$ (partial^3 H)/(partial x^3)|_(x = L) = 0 => 1/2 H(L + 2 Delta x) - H(L + Delta x) + H(L - Delta x) - 1/2 H(L - 2 Delta x) = 0 $
$ => H(L + 2 Delta x) = 2 H(L + Delta x) - 2 H(L - Delta x) + H(L - 2 Delta x) $

== Frequencies

Suppose $H(x, t) = h(x) cos(omega t)$. Plugging into the beam equation:
$ E I (d^4 h)/(d x^4) cos(omega t) = rho A h(x) omega^2 cos(omega t) quad "[ODE]" $

Solution:
$ f_n = omega_n / (2 pi) = (beta_n L)^2 / (2 pi L^2) sqrt(E I / (rho A)) quad "+ boundary conditions" $

where $beta_1 L approx 1.875104$.

== CFL Conditions

$ Delta t <= (Delta x^2) / 2 sqrt(rho A / (E I)) $

== Damping

$ cases(
  (partial U)/(partial t) = -(E I)/(rho A) (partial^4 H)/(partial x^4) - gamma U quad "(damping)",
  (partial H)/(partial t) = U quad "(unchanged)"
) $

= Results

== Simulation Parameters

$ E_"iron" = 2.1 times 10^11 "Pa", quad rho_"iron" = 7874 "kg/m"^3 $

$ F = 3 "N", quad gamma = 8 $

Beam radius and cross-section:
$ R = 0.5 times 10^(-3) "m", quad A = pi R^2, quad I = pi R^4 / 4 $

Frequencies:
$ f_n = 256 dot 2^(n / 12), quad n = 0, 1, ..., 11 $

Audio: #emph[[audio file]]

= References

- https://ocw.mit.edu/courses/2-002-mechanics-and-materials-ii-spring-2004/bc25a56b5a91ad29ca5c7419616686f7_lec2.pdf
- https://web.ics.purdue.edu/~gonza226/ME323/Lecture-16.pdf
- Hibbeler, _Mechanics of Materials_ (10th Edition)
