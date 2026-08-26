//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "RadialReturnStressUpdate.h"

/**
 * Base class for isotropic scalar radial-return models whose residual has the form of a
 * rate-dependent inelastic strain increment minus the return-mapping scalar.
 */
template <bool is_ad>
class RadialReturnRateDependentStressUpdateTempl : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  RadialReturnRateDependentStressUpdateTempl(const InputParameters & parameters);

protected:
  enum class InitialGuessType
  {
    ZERO,
    WEN,
    OLD_STRESS
  };

  virtual GenericReal<is_ad>
  initialGuess(const GenericReal<is_ad> & effective_trial_stress) override;

  /**
   * Compute the Equation 13 initial guess using the model state at the old stress.
   * @param effective_trial_stress Effective trial stress
   * @param wen_guess Scalar that reduces the trial stress to the old von Mises stress
   */
  virtual GenericReal<is_ad>
  oldStressInitialGuess(const GenericReal<is_ad> & effective_trial_stress,
                        const GenericReal<is_ad> & wen_guess);

  /// Method used to initialize the non-AD scalar return-mapping solve
  const InitialGuessType _initial_guess_type;
};

typedef RadialReturnRateDependentStressUpdateTempl<false> RadialReturnRateDependentStressUpdate;
typedef RadialReturnRateDependentStressUpdateTempl<true> ADRadialReturnRateDependentStressUpdate;
