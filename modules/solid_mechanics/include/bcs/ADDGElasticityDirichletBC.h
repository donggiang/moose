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
 * Weak Dirichlet boundary condition for vector linear elasticity in
 * an SIPG / NIPG / IIPG / OBB DG formulation, following Liu, Wheeler,
 * Dawson (2009) Eq. 22 boundary contribution and Rabedi note 8 eq. 6
 * (Dirichlet specialisation: u^* = u_bar at Dirichlet face).
 *
 * Boundary face residual for component k of the displacement:
 *   r = - (sigma(u) . n)_k                                (consistency)
 *       + epsilon * (sigma(v) . n)_j (u - u_bar)_j        (adjoint)
 *       + (eta * G / h_F) (u - u_bar)_k                   (penalty)
 *
 * with v = e_k * test_scalar (only kernel's component non-zero), and
 * u_bar provided per-component by a Function.
 *
 * One BC instance per displacement component, like the existing
 * solid-mechanics BC pattern. Each instance reads a separate Function
 * for the prescribed displacement of that component.
 */
class ADDGElasticityDirichletBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ADDGElasticityDirichletBC(const InputParameters & parameters);

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

  /// Component (0/1/2) the BC acts on
  const unsigned int _component;
  const unsigned int _ndisp;

  /// Penalty parameters (Liu 2009)
  const Real _eta;
  const Real _epsilon;

  /// Per-component prescribed displacement functions
  std::vector<const Function *> _u_bar;

  /// Element materials
  const ADMaterialProperty<RankTwoTensor> & _stress;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  const ADMaterialProperty<Real> & _penalty_modulus;

  /// All displacement variables, for vector deviation (u - u_bar)
  std::vector<const ADVariableValue *> _disp;
};
