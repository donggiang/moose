//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideIntegralPostprocessor.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"

class Function;

/**
 * Discrete reaction extractor for the Abedi star-value DG Dirichlet BC.
 *
 * Mirrors the bilinear form contribution that ADDGElasticityAbediDirichletBC
 * adds on a Dirichlet sideset, but with the FE test function replaced by a
 * user-prescribed virtual displacement field v(x) = direction * f(x):
 *
 *   R = INT_Gamma  [ - (sigma(u_h) . n) . v
 *                    + epsilon * (sigma(v) . n) . (u_h - u_bar) ] dS
 *
 *   - direction (required): unit vector picking out the reaction component
 *   - f(x) (`test_lifting`, optional): scalar lifting; default f = 1
 *   - epsilon: DG variant (must match the BC's epsilon)
 *   - function: u_bar (per displacement component); default 0
 *
 * For f = 1 (constant) the adjoint term vanishes (sigma(v)=0) and the
 * integral collapses to the surface stress integral -- equivalent to
 * ADSidesetReaction (up to sign). The lifting parameter exists to make the
 * adjoint contribution explicit and inspectable.
 */
class ADDGAbediSidesetReaction : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();

  ADDGAbediSidesetReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const unsigned int _ndisp;
  const RealVectorValue _direction;
  const Real _epsilon;

  const ADMaterialProperty<RankTwoTensor> & _stress;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;

  std::vector<const ADVariableValue *> _disp;
  std::vector<const Function *> _u_bar;
  const Function * const _test_lifting;
};
