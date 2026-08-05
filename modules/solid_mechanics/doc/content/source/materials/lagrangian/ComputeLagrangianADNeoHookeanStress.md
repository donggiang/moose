# ComputeLagrangianADNeoHookeanStress

!syntax description /Materials/ComputeLagrangianADNeoHookeanStress

## Overview

`ComputeLagrangianADNeoHookeanStress` demonstrates selective automatic differentiation within
the non-AD Lagrangian mechanics workflow. The existing Lagrangian implementation computes the
Green-Lagrange strain, Neo-Hookean PK2 stress, stress-measure transformations, and kinematic
Jacobian terms. Local forward automatic differentiation is used only to compute the material
tangent $\partial S_{ij}/\partial E_{kl}$ from the PK2 stress relation.

The ordinary PK2 stress is

\begin{equation}
  \boldsymbol{S} = (\lambda \log J - \mu)\boldsymbol{C}^{-1} + \mu\boldsymbol{I},
\end{equation}

where $\boldsymbol{C}=2\boldsymbol{E}+\boldsymbol{I}$. A local AD copy of this same relation
seeds the nine components of $\boldsymbol{E}$ and extracts

\begin{equation}
  C_{ijkl} = \frac{\partial S_{ij}}{\partial E_{kl}}.
\end{equation}

The local derivatives are stripped into an ordinary `RankFourTensor`. The existing non-AD
`ComputeLagrangianStressPK2` base then converts PK2 to PK1 and chains the material tangent with
the existing kinematic derivatives before the non-AD kernel and assembly system consume it.

!alert note
This material supports large kinematics only.

## Example Input File Syntax

!listing modules/solid_mechanics/test/tests/lagrangian/materials/convergence/neohookean.i
         block=Materials

!syntax parameters /Materials/ComputeLagrangianADNeoHookeanStress

!syntax inputs /Materials/ComputeLagrangianADNeoHookeanStress

!syntax children /Materials/ComputeLagrangianADNeoHookeanStress
