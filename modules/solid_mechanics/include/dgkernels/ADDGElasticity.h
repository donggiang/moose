//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDGKernel.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"

/**
 * SIPG/NIPG/IIPG/OBB DG kernel for vector linear elasticity, following
 * Liu, Wheeler, Dawson (2009) Eq. 22 and Rabedi note 7. Acts on a single
 * displacement component; one instance per component is required for
 * full vector elasticity.
 *
 * Interior face residual contribution:
 *   r = - <sigma(u) . n>_i * [v]_i              (consistency)
 *       + epsilon * <sigma(v) . n>_j * [u]_j    (adjoint)
 *       + (eta * G / h_F) * [u]_i * [v]_i       (penalty, Liu form)
 *
 * where
 *   epsilon  = -1 SIPG, +1 NIPG, 0 IIPG (and OBB has eta = 0)
 *   eta      = dimensionless penalty parameter (Liu's delta_p)
 *   G        = penalty modulus (typically shear modulus mu, but read
 *              from a material property so the user picks the
 *              dimensional scale)
 *   h_F      = face characteristic size
 *
 * The adjoint term involves sigma(v) for a hypothetical test
 * displacement v = e_i v_i (only the kernel's component is nonzero).
 * It is computed via the elasticity tensor C: sigma(v) = C : eps(v).
 * This term couples the test function (one component) to the FULL
 * jump in displacement [u]_j across all components.
 */
class ADDGElasticity : public ADDGKernel
{
public:
  static InputParameters validParams();
  ADDGElasticity(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  /// Helper: build full rank-2 displacement-gradient tensor for the
  /// scalar test function in component _component (the only non-zero row).
  /// Templated on input type to accept both Real and AD vectors.
  template <typename T>
  ADRankTwoTensor buildGradV(const T & grad_test_scalar) const
  {
    ADRankTwoTensor grad_v;
    for (unsigned int p = 0; p < LIBMESH_DIM; ++p)
      grad_v(_component, p) = grad_test_scalar(p);
    return grad_v;
  }

  /// Component (0/1/2) the kernel acts on
  const unsigned int _component;

  /// Number of displacement components in the problem
  const unsigned int _ndisp;

  /// Dimensionless penalty parameter (Liu 2009 delta_p)
  const Real _eta;

  /// DG variant indicator: -1 SIPG, 0 IIPG, +1 NIPG
  const Real _epsilon;

  /// Stress tensor on this side and neighbor (consistency)
  const ADMaterialProperty<RankTwoTensor> & _stress;
  const ADMaterialProperty<RankTwoTensor> & _stress_neighbor;

  /// Elasticity tensor on this side and neighbor (adjoint: sigma(v) = C:eps(v))
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor_neighbor;

  /// Penalty modulus (Liu's G); name configurable via input parameter
  const ADMaterialProperty<Real> & _penalty_modulus;
  const ADMaterialProperty<Real> & _penalty_modulus_neighbor;

  /// All displacement components, this side (for vector jump [u]_j)
  std::vector<const ADVariableValue *> _disp;
  std::vector<const ADVariableValue *> _disp_neighbor;
};
