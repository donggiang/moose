# ComputeLagrangianADNeoHookeanStress

!syntax description /Materials/ComputeLagrangianADNeoHookeanStress

## Overview

`ComputeLagrangianADNeoHookeanStress` demonstrates selective automatic differentiation within
the non-AD Lagrangian mechanics workflow. Nested local forward automatic differentiation seeds
the nine components of the deformation gradient, computes the right Cauchy-Green tensor
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
material properties before the existing non-AD kernel and assembly system consume them. The AD
types therefore remain confined to the material calculation.

!alert note
This material supports large kinematics only.

## Example Input File Syntax

!listing modules/solid_mechanics/test/tests/lagrangian/materials/convergence/neohookean.i
         block=Materials

!syntax parameters /Materials/ComputeLagrangianADNeoHookeanStress

!syntax inputs /Materials/ComputeLagrangianADNeoHookeanStress

!syntax children /Materials/ComputeLagrangianADNeoHookeanStress
