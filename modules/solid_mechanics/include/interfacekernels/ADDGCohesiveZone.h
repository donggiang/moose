//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADInterfaceKernel.h"

/**
 * DG cohesive-zone interface kernel using the Abedi/Rabedi star-value
 * pattern with a Neumann-type substitution on the cohesive face:
 *   t*.n_e = T_coh([[u]]) . n_e     (TSL replaces the trial flux)
 *   u* = u                          (no replacement; adjoint vanishes)
 *
 * Bilinear contribution at a cohesive face Gamma_c, summed over both
 * sides:
 *   - INT_Gamma_c [[w]] . T_coh dS
 *
 * In MOOSE InterfaceKernel form (residual per side):
 *   r_element  = -T_coh_on_element[c] * test
 *   r_neighbor = +T_coh_on_element[c] * test_neighbor
 *
 * The constitutive law is a bilinear traction-separation with
 * Benzeggagh-Kenane (B-K) mixed-mode interpolation (Turon 2006). The
 * "minus" branch on the normal direction (closing) is kept elastic --
 * compression is not damaged. Damage is computed from the CURRENT
 * effective separation; this implementation does NOT track history,
 * which is acceptable for monotonic loading (e.g. mode-I DCB).
 *
 * One instance per displacement component is required.
 */
class ADDGCohesiveZone : public ADInterfaceKernel
{
public:
  static InputParameters validParams();
  ADDGCohesiveZone(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  /// Displacement component this kernel acts on (0=x, 1=y, 2=z)
  const unsigned int _component;

  /// Number of displacement components in the problem
  const unsigned int _ndisp;

  /// Cohesive parameters (Turon-style bilinear with B-K)
  const Real _K;       ///< penalty stiffness [F/L^3]
  const Real _T;       ///< peak normal traction [F/L^2]
  const Real _S;       ///< peak shear traction [F/L^2]
  const Real _GIc;     ///< mode-I fracture energy [F/L]
  const Real _GIIc;    ///< mode-II fracture energy [F/L]
  const Real _eta_BK;  ///< Benzeggagh-Kenane mixed-mode parameter

  /// Optional viscous regularization (Kelvin-Voigt damper in parallel with
  /// the cohesive TSL). T_visc = -viscosity * d/dt [[u]]. Smooths the
  /// post-peak softening branch so Newton can traverse the snap. Set 0 to
  /// recover the rate-independent law.
  const Real _viscosity;

  /// Residual stiffness floor as a fraction of K, used to avoid singular
  /// tangent at full damage. Effective degradation factor becomes
  /// (1-d)*(1-kappa) + kappa, so at d=1 the tangent stays kappa*K instead
  /// of dropping to zero. Set 0 for the classical fully-degraded TSL.
  const Real _kappa;

  /// All displacement components on each side (for vector jump)
  std::vector<const ADVariableValue *> _disp;
  std::vector<const ADVariableValue *> _disp_neighbor;

  /// Old (previous converged step) displacement values on each side for
  /// the viscous rate term. These are non-AD; viscous Jacobian comes from
  /// the current jump via AD.
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableValue *> _disp_neighbor_old;
};
