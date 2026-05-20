# Tier 1.1 — Mathematical formulation

This document records the equations underlying
`tier1_1_mode_I_monotonic.i` and the analytical reference computed in
`verify_tier1_1.py`.

## Domain and discontinuity

Two stacked unit-square subdomains $\Omega_1=[0,1]\times[0,1]$ and
$\Omega_2=[0,1]\times[1,2]$, with `BreakMeshByBlockGenerator` duplicating the
nodes along the shared edge $\Gamma_c=\{(x,y):0\le x\le 1,\,y=1\}$. The
displacement field is allowed to be discontinuous across $\Gamma_c$:

$$
\llbracket \mathbf{u} \rrbracket
= \mathbf{u}^{+} - \mathbf{u}^{-}\quad \text{on } \Gamma_c,
$$

where the superscripts denote the two sides of $\Gamma_c$.
The interface unit normal $\mathbf{n}=(0,1)$ points from $\Omega_1$ into
$\Omega_2$. The normal and tangential components of the jump are

$$
\delta_n = \llbracket \mathbf{u} \rrbracket \cdot \mathbf{n}
= \llbracket u_y \rrbracket, \qquad
\delta_t = \llbracket \mathbf{u} \rrbracket \cdot \mathbf{t}
= \llbracket u_x \rrbracket .
$$

## Bulk equations (small strain, linear elastic)

In each subdomain we solve the static momentum balance with no body force,

$$
\nabla\!\cdot \boldsymbol{\sigma} = \mathbf{0}\quad\text{in }\Omega_i,\;i=1,2,
$$

with isotropic linear elasticity

$$
\boldsymbol{\sigma}
= \mathbb{C}:\boldsymbol{\varepsilon},\qquad
\boldsymbol{\varepsilon}
= \tfrac{1}{2}\bigl(\nabla\mathbf{u}+\nabla\mathbf{u}^{\!\top}\bigr),
$$

and isotropic constants $E=10^{9}$, $\nu=0.3$. The bulk is intentionally
much stiffer than the interface so the global displacement is dominated by
the cohesive jump.

## Interface kinematics, traction and balance

On $\Gamma_c$ the cohesive interface kernel enforces

$$
\boldsymbol{\sigma}\,\mathbf{n} = \mathbf{t}(\delta_n,\delta_t),
$$

where $\mathbf{t}$ is the traction supplied by the cohesive material. In
mode I the only active component is

$$
t_n = \mathbf{t}\cdot\mathbf{n}.
$$

## Bilinear (Camanho–Dávila) traction–separation law

In Tier 1.1 the loading is monotonic and pure mode I, so $\delta_t=0$ and
the maximum jump equals the current jump,
$\delta_n^{\max}=\delta_n$. The bilinear law of
[Camanho & Dávila, NASA/TM-2002-211737] reduces to:

$$
t_n(\delta_n) =
\begin{cases}
K\,\delta_n
& 0 \le \delta_n \le \delta_n^{0}, \\[6pt]
(1-d)\,K\,\delta_n
& \delta_n^{0} < \delta_n \le \delta_n^{f}, \\[6pt]
0
& \delta_n > \delta_n^{f},
\end{cases}
$$

with the **scalar damage variable**

$$
d \;=\; \frac{\delta_n^{f}\,\bigl(\delta_n^{\max} - \delta_n^{0}\bigr)}
              {\delta_n^{\max}\,\bigl(\delta_n^{f} - \delta_n^{0}\bigr)},
\qquad d\in[0,1].
$$

For monotonic loading the softening branch can be written equivalently in
the linear form used by `verify_tier1_1.py`:

$$
t_n(\delta_n) = N\,\frac{\delta_n^{f} - \delta_n}{\delta_n^{f} - \delta_n^{0}},
\quad \delta_n^{0} \le \delta_n \le \delta_n^{f}.
$$

The two characteristic separations follow from the peak-strength and
fracture-energy conditions:

$$
\delta_n^{0} = \frac{N}{K}, \qquad
\delta_n^{f} = \frac{2\,G_{Ic}}{N}.
$$

The latter is the area-under-the-curve identity

$$
G_{Ic} \;=\; \int_{0}^{\delta_n^{f}} t_n(\delta_n)\,d\delta_n
\;=\; \tfrac{1}{2}\,N\,\delta_n^{f}.
$$

## Numerical values used in Tier 1.1

| Symbol | Meaning | Value |
|---|---|---|
| $K$ | penalty stiffness | $1.0\times 10^{6}$ |
| $N$ | normal strength | $50$ |
| $G_{Ic}$ | mode-I fracture energy | $0.5$ |
| $\delta_n^{0}=N/K$ | onset of softening | $5.0\times 10^{-5}$ |
| $\delta_n^{f}=2G_{Ic}/N$ | full decohesion | $2.0\times 10^{-2}$ |

## Loading and boundary conditions

The bottom face $y=0$ is fully clamped; the top face $y=2$ is constrained
in $x$ and pulled in $y$ by a piecewise-linear function $u_*(t)$ that
spends 10 % of the simulation time in the elastic regime so both branches
are well sampled by a uniform $\Delta t=0.01$:

$$
u_*(t) = \begin{cases}
\dfrac{7.5\times 10^{-5}}{0.1}\;t,
& 0 \le t \le 0.1, \\[8pt]
7.5\times 10^{-5}
+ \dfrac{0.03 - 7.5\times 10^{-5}}{0.9}\,(t-0.1),
& 0.1 < t \le 1.0 .
\end{cases}
$$

The endpoint $u_*(1)=0.03=1.5\,\delta_n^{f}$ guarantees we run past full
decohesion.

## Discrete weak form (sketch)

Find $\mathbf{u}\in\mathcal{V}$ such that for all
$\mathbf{v}\in\mathcal{V}_0$,

$$
\underbrace{\int_{\Omega_1\cup\Omega_2}\!
\boldsymbol{\sigma}(\mathbf{u}):\nabla \mathbf{v}\;dV}_{\text{bulk}}
\;+\;
\underbrace{\int_{\Gamma_c}
\mathbf{t}(\llbracket\mathbf{u}\rrbracket)\cdot
\llbracket\mathbf{v}\rrbracket\;dS}_{\text{cohesive interface}}
\;=\;0 .
$$

The cohesive integral is implemented in MOOSE as `CZMInterfaceKernelSmallStrain`
(activated by `Physics/SolidMechanics/CohesiveZone` with `strain = SMALL`).
At every interface quadrature point, $\mathbf{t}$ is the value returned by
`BiLinearMixedModeTraction::computeTraction()` and its consistent tangent
$\partial \mathbf{t}/\partial \llbracket \mathbf{u}\rrbracket$ is the
hand-coded `computeTractionDerivatives()`.

## Pass criterion

Let $\{(\delta_n^{(k)},\, t_n^{(k)})\}_{k=0}^{100}$ be the simulated
samples and $t_n^{\text{ref}}(\delta)$ the closed-form law above. The
verification accepts the run if

$$
\max_{k}\;\bigl|\,t_n^{(k)} - t_n^{\text{ref}}(\delta_n^{(k)})\,\bigr|
\;<\;10^{-7}, \qquad
\max_{k}\;|t_t^{(k)}|\;<\;10^{-10}.
$$

The first condition checks the mode-I curve; the second confirms that no
spurious tangential traction develops under pure normal loading. For the
parameters above the simulation gives

$$
\max_k\;\bigl|t_n - t_n^{\text{ref}}\bigr| = 4.5\times 10^{-12},\qquad
\max_k\;|t_t| = 3.3\times 10^{-14},
$$

i.e. the law is reproduced to roughly 13 significant figures relative to
the peak traction $N=50$.
