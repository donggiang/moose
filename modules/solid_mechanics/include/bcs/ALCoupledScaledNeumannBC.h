//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

class Function;

/**
 * Neumann BC scaled by a coupled scalar variable λ:
 *
 *     residual contribution = − λ · f_ext(x, t) · φ_test
 *
 * where `f_ext` comes from a `function` parameter (same convention as
 * `FunctionNeumannBC`) and `λ` is the scalar variable provided by
 * `arc_length_scalar`. Off-diagonal Jacobian wrt λ is filled so that the
 * augmented arc-length system has a complete bordered tangent.
 *
 * Companion to `ArcLengthScalarKernel`. Together they implement Path C
 * (input-only augmented system) arc-length: PETSc's standard Newton solves
 * the bordered system without needing SNESNEWTONAL.
 */
class ALCoupledScaledNeumannBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  ALCoupledScaledNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override { return 0.0; }
  virtual Real computeQpOffDiagJacobianScalar(unsigned int svar) override;

  /// External-load function (the q^ in arc-length parlance).
  const Function & _func;

  /// Coupled scalar variable values (λ); index [0] for the scalar's value.
  const VariableValue & _lambda;

  /// Scalar variable's variable number (to discriminate in off-diag Jacobian).
  const unsigned int _lambda_var;
};
