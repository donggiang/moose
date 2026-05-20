//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KernelScalarBase.h"

/**
 * Cylindrical (Riks) arc-length constraint kernel:
 *
 *     g(u, λ) = ∫_Ω |Δu|^2 dΩ  -  Δs^2  =  0
 *
 * where Δu = u - u_old (component-wise), λ is the load parameter
 * (a scalar variable kept here as `_kappa`), and Δs is the prescribed
 * arc-length increment per timestep.
 *
 * One instance of this kernel is required per displacement component.
 * The component's per-QP contribution to the constraint is `(_u - _u_old)^2`.
 * Exactly one instance must be designated `primary_constraint = true`,
 * which additionally subtracts `Δs^2 / V_domain` per QP so that the integral
 * yields the global `-Δs^2` term.
 *
 * Inputs needed:
 *   - variable        : the displacement component this kernel focuses on
 *   - scalar_variable : the load-parameter scalar variable (commonly `lambda`)
 *   - delta_s         : prescribed arc-length step (Δs)
 *   - volume_pp       : Postprocessor giving V_domain (e.g. VolumePostprocessor)
 *   - primary_constraint : true on exactly ONE instance per λ; false on the others
 *
 * The off-diagonal Jacobian d-(constraint)/d-(displacement) is filled with
 * 2 (u - u_old) ∂u/∂u_j = 2 (_u - _u_old) _phi[_j][_qp]. The diagonal Jacobian
 * d-(constraint)/d-λ is zero in the cylindrical formulation.
 */
class ArcLengthScalarKernel : public KernelScalarBase
{
public:
  static InputParameters validParams();
  ArcLengthScalarKernel(const InputParameters & parameters);

protected:
  /// Field-residual contribution from this kernel: zero (the load BC handles -λ·F_ext).
  virtual Real computeQpResidual() override { return 0.0; }
  virtual Real computeQpJacobian() override { return 0.0; }

  /// Scalar (constraint) residual integrand at this QP.
  virtual Real computeScalarQpResidual() override;

  /// Off-diagonal Jacobian d-(constraint)/d-(field var) at this QP.
  virtual Real computeScalarQpOffDiagJacobian(const unsigned int jvar_num) override;

  /// Diagonal Jacobian d-(constraint)/d-λ — zero for cylindrical formulation.
  virtual Real computeScalarQpJacobian() override { return 0.0; }

  /// Prescribed arc-length step.
  const Real _delta_s;

  /// Total domain volume (for distributing the -Δs^2 constant inside the QP integrand).
  const PostprocessorValue & _volume;

  /// True if this instance carries the -Δs^2 constant.
  const bool _primary_constraint;

  /// Old solution at this QP (component-wise).
  const VariableValue & _u_old_qp;
};
