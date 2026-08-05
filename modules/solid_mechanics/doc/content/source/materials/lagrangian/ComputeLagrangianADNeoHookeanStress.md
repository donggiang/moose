# ComputeLagrangianADNeoHookeanStress

!syntax description /Materials/ComputeLagrangianADNeoHookeanStress

## Overview

`ComputeLagrangianADNeoHookeanStress` demonstrates selective automatic differentiation within
the non-AD Lagrangian mechanics workflow. Nested local forward automatic differentiation starts
from the displacement gradient and constructs

\begin{equation}
  \boldsymbol{F}=\boldsymbol{I}+\nabla_X\boldsymbol{u}.
\end{equation}

It then computes the right Cauchy-Green tensor
$\boldsymbol{C}=\boldsymbol{F}^T\boldsymbol{F}$, and evaluates the compressible Neo-Hookean
strain-energy potential

\begin{equation}
  \Psi = \frac{\lambda}{2}(\log J)^2 - \mu\log J
       + \frac{\mu}{2}(\operatorname{tr}\boldsymbol{C}-3),
\end{equation}

Its gradient with respect to the deformation gradient gives the PK1 stress,

\begin{equation}
  P_{ij} = \frac{\partial \Psi}{\partial F_{ij}},
\end{equation}

and its Hessian includes both constitutive and kinematic derivatives in the material tangent,

\begin{equation}
  \mathcal{A}_{ijkl} = \frac{\partial^2 \Psi}{\partial F_{ij}\partial F_{kl}}.
\end{equation}

The local derivatives are stripped into ordinary `RankTwoTensor` and `RankFourTensor` PK1
material properties before the existing non-AD kernel consumes them. The kernel retains the
explicit residual and stiffness expressions: it contracts stress with the test-function gradient
and the material tangent with the test- and trial-function gradients. Thus the final chain from
$\nabla_X\boldsymbol{u}=\sum_b\boldsymbol{u}_b\otimes\nabla_X N_b$ to nodal displacement is
assembled explicitly, while all constitutive and deformation-gradient operations remain in the
local AD graph.

For Cartesian total-Lagrangian mechanics, the resulting element terms are

\begin{equation}
  R_a = \int_{\Omega_0} \nabla_X N_a : \boldsymbol{P}\,dV,
\end{equation}

and

\begin{equation}
  K_{ab} = \int_{\Omega_0} \nabla_X N_a : \mathcal{A} : \nabla_X N_b\,dV,
  \qquad \mathcal{A}=\frac{\partial\boldsymbol{P}}{\partial\boldsymbol{F}}.
\end{equation}

These are assembled by `TotalLagrangianStressDivergence`; this material does not replace that
non-AD kernel.

!alert note
This material supports large kinematics only.

## Example Input File Syntax

!listing modules/solid_mechanics/test/tests/lagrangian/materials/convergence/neohookean.i
         block=Materials

!syntax parameters /Materials/ComputeLagrangianADNeoHookeanStress

!syntax inputs /Materials/ComputeLagrangianADNeoHookeanStress

!syntax children /Materials/ComputeLagrangianADNeoHookeanStress
