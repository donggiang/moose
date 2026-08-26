//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "RadialReturnRateDependentStressUpdate.h"

#include "MathUtils.h"
#include "MooseUtils.h"

template <bool is_ad>
InputParameters
RadialReturnRateDependentStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  MooseEnum initial_guess_type("ZERO WEN OLD_STRESS", "OLD_STRESS");
  params.addParam<MooseEnum>(
      "initial_guess_type",
      initial_guess_type,
      "Initial guess for the scalar return-mapping solve. ZERO starts from a zero effective "
      "inelastic strain increment; WEN estimates the increment from the difference between the "
      "trial and old von Mises stresses; OLD_STRESS (the default) evaluates the rate-dependent "
      "inelastic increment using the old stress and internal state.");
  return params;
}

template <bool is_ad>
RadialReturnRateDependentStressUpdateTempl<is_ad>::
    RadialReturnRateDependentStressUpdateTempl(const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),
    _initial_guess_type(
        parameters.get<MooseEnum>("initial_guess_type").template getEnum<InitialGuessType>())
{
}

template <bool is_ad>
GenericReal<is_ad>
RadialReturnRateDependentStressUpdateTempl<is_ad>::initialGuess(
    const GenericReal<is_ad> & effective_trial_stress)
{
  if (_initial_guess_type == InitialGuessType::ZERO ||
      MooseUtils::absoluteFuzzyEqual(this->_effective_old_stress, 0.0))
    return 0.0;

  const GenericReal<is_ad> wen_guess =
      (effective_trial_stress - this->_effective_old_stress) / this->_three_shear_modulus;
  const GenericReal<is_ad> initial_guess =
      _initial_guess_type == InitialGuessType::WEN
          ? wen_guess
          : oldStressInitialGuess(effective_trial_stress, wen_guess);

  return MathUtils::clamp(initial_guess,
                          this->minimumPermissibleValue(effective_trial_stress),
                          this->maximumPermissibleValue(effective_trial_stress));
}

template <bool is_ad>
GenericReal<is_ad>
RadialReturnRateDependentStressUpdateTempl<is_ad>::oldStressInitialGuess(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & wen_guess)
{
  return wen_guess + this->computeResidual(effective_trial_stress, wen_guess);
}

template class RadialReturnRateDependentStressUpdateTempl<false>;
template class RadialReturnRateDependentStressUpdateTempl<true>;
