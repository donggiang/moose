# ComputeLagrangianADNeoHookeanStress

!syntax description /Materials/ComputeLagrangianADNeoHookeanStress

## Overview

`ComputeLagrangianADNeoHookeanStress` demonstrates selective automatic differentiation within
the non-AD Lagrangian mechanics workflow. The existing Lagrangian implementation computes the
Green-Lagrange strain, stress-measure transformations, and kinematic Jacobian terms. Nested local
forward automatic differentiation computes the PK2 stress and material tangent from the
compressible Neo-Hookean strain-energy potential

\begin{equation}
  \Psi = \frac{\lambda}{2}(\log J)^2 - \mu\log J
       + \frac{\mu}{2}(\operatorname{tr}\boldsymbol{C}-3),
\end{equation}

where $\boldsymbol{C}=2\boldsymbol{E}+\boldsymbol{I}$. Its gradient gives the PK2 stress,

\begin{equation}
  S_{ij} = \frac{\partial \Psi}{\partial E_{ij}},
\end{equation}

and its Hessian gives the material tangent,

\begin{equation}
  C_{ijkl} = \frac{\partial^2 \Psi}{\partial E_{ij}\partial E_{kl}}.
\end{equation}

The local derivatives are stripped into ordinary `RankTwoTensor` and `RankFourTensor` material
properties. The existing non-AD `ComputeLagrangianStressPK2` base then converts PK2 to PK1 and
chains the material tangent with the existing kinematic derivatives before the non-AD kernel and
assembly system consume them.

!alert note
This material supports large kinematics only.

## Example Input File Syntax

!listing modules/solid_mechanics/test/tests/lagrangian/materials/convergence/neohookean.i
         block=Materials

!syntax parameters /Materials/ComputeLagrangianADNeoHookeanStress

!syntax inputs /Materials/ComputeLagrangianADNeoHookeanStress

!syntax children /Materials/ComputeLagrangianADNeoHookeanStress
