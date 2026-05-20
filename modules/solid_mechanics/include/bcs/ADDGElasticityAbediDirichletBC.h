//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"

class Function;

/**
 * Abedi star-value Dirichlet BC for vector linear elasticity in DG
 * (elasticity analogue of DGAbediDirichletBC for the heat equation).
 * Substitutes u* = u_bar and (sigma.n)* = sigma.n (no replacement)
 * into the general DG bilinear form (Abedi/Rabedi note eq. 113/115),
 * keeping only consistency + adjoint-consistency on Dirichlet faces.
 *
 * NO G/h_F penalty term -- that piece is what ADDGElasticityDirichletBC
 * adds for SIPG coercivity but is OPTIONAL in the abstract derivation.
 *
 * Boundary face residual for component k of the displacement:
 *   r = - (sigma(u) . n)_k * test                       (consistency)
 *       + epsilon * (sigma(v) . n) . (u - u_bar)        (adjoint)
 *
 * with v = e_k * test_scalar (only the BC's component non-zero), and
 * u_bar provided per-component by a Function.
 *
 * One BC instance per displacement component. NIPG (epsilon = +1) is
 * stable without the penalty; SIPG (epsilon = -1) without the penalty
 * is non-coercive and will produce a divergent/garbage solution.
 */
class ADDGElasticityAbediDirichletBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ADDGElasticityAbediDirichletBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Build full rank-2 grad-of-test tensor for scalar test on _component.
  template <typename T>
  ADRankTwoTensor buildGradV(const T & grad_test_scalar) const
  {
    ADRankTwoTensor grad_v;
    for (unsigned int p = 0; p < LIBMESH_DIM; ++p)
      grad_v(_component, p) = grad_test_scalar(p);
    return grad_v;
  }

  const unsigned int _component;
  const unsigned int _ndisp;
  const Real _epsilon;

  std::vector<const Function *> _u_bar;

  const ADMaterialProperty<RankTwoTensor> & _stress;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;

  std::vector<const ADVariableValue *> _disp;
};
