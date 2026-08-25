//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PowerLawCreepStressUpdate.h"

#include "MathUtils.h"
#include "MooseUtils.h"

registerMooseObject("SolidMechanicsApp", PowerLawCreepStressUpdate);
registerMooseObject("SolidMechanicsApp", ADPowerLawCreepStressUpdate);

template <bool is_ad>
InputParameters
PowerLawCreepStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic power law creep "
      "model. This class can be used in conjunction with other creep and plasticity materials "
      "for more complex simulations.");

  // Linear strain hardening parameters
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addRequiredParam<Real>("coefficient", "Leading coefficient in power-law equation");
  params.addRequiredParam<Real>("n_exponent", "Exponent on effective stress in power-law equation");
  params.addParam<Real>("m_exponent", 0.0, "Exponent on time in power-law equation");
  params.addRequiredParam<Real>("activation_energy", "Activation energy");
  params.addParam<Real>("gas_constant", 8.3143, "Universal gas constant");
  params.addParam<Real>("start_time", 0.0, "Start time (if not zero)");
  MooseEnum initial_guess_type("ZERO WEN OLD_STRESS", "WEN");
  params.addParam<MooseEnum>(
      "initial_guess_type",
      initial_guess_type,
      "Initial guess for the non-AD scalar return-mapping solve. ZERO starts from a zero "
      "effective creep strain increment; WEN (the default) estimates the increment from the "
      "difference between the trial and old von Mises stresses; OLD_STRESS evaluates the "
      "power-law creep increment "
      "using the old von Mises stress. AD solves retain the zero initial guess.");
  return params;
}

template <bool is_ad>
PowerLawCreepStressUpdateTempl<is_ad>::PowerLawCreepStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _initial_guess_type(
        parameters.get<MooseEnum>("initial_guess_type").template getEnum<InitialGuessType>()),
    _temperature(this->isParamValid("temperature")
                     ? &this->template coupledGenericValue<is_ad>("temperature")
                     : nullptr),
    _coefficient(this->template getParam<Real>("coefficient")),
    _n_exponent(this->template getParam<Real>("n_exponent")),
    _m_exponent(this->template getParam<Real>("m_exponent")),
    _activation_energy(this->template getParam<Real>("activation_energy")),
    _gas_constant(this->template getParam<Real>("gas_constant")),
    _start_time(this->template getParam<Real>("start_time")),
    _exponential(1.0)
{
  if (_start_time < this->_app.getStartTime() && (std::trunc(_m_exponent) != _m_exponent))
    this->paramError("start_time",
                     "Start time must be equal to or greater than the Executioner start_time if a "
                     "non-integer m_exponent is used");
}

template <bool is_ad>
void
PowerLawCreepStressUpdateTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  using std::exp, std::pow;

  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  if (_temperature)
    _exponential = exp(-_activation_energy / (_gas_constant * (*_temperature)[_qp]));

  _exp_time = pow(_t - _start_time, _m_exponent);
}

template <bool is_ad>
GenericReal<is_ad>
PowerLawCreepStressUpdateTempl<is_ad>::initialGuess(
    const GenericReal<is_ad> & effective_trial_stress)
{
  // An AD initial guess can be accepted before a Newton update establishes the derivative of the
  // converged constitutive response, so retain the established zero guess for AD solves.
  if constexpr (is_ad)
    return 0.0;

  if (_initial_guess_type == InitialGuessType::ZERO ||
      MooseUtils::absoluteFuzzyEqual(this->_effective_old_stress, 0.0))
    return 0.0;

  GenericReal<is_ad> initial_guess;
  if (_initial_guess_type == InitialGuessType::WEN)
    initial_guess =
        (effective_trial_stress - this->_effective_old_stress) / this->_three_shear_modulus;
  else
  {
    using std::pow;
    initial_guess = _dt * _coefficient * pow(this->_effective_old_stress, _n_exponent) *
                    _exponential * _exp_time;
  }

  return MathUtils::clamp(initial_guess,
                          this->minimumPermissibleValue(effective_trial_stress),
                          this->maximumPermissibleValue(effective_trial_stress));
}

template <bool is_ad>
template <typename ScalarType>
ScalarType
PowerLawCreepStressUpdateTempl<is_ad>::computeResidualInternal(
    const GenericReal<is_ad> & effective_trial_stress, const ScalarType & scalar)
{
  using std::pow;

  const ScalarType stress_delta = effective_trial_stress - _three_shear_modulus * scalar;
  const ScalarType creep_rate =
      _coefficient * pow(stress_delta, _n_exponent) * _exponential * _exp_time;
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
PowerLawCreepStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  using std::pow;

  const GenericReal<is_ad> stress_delta = effective_trial_stress - _three_shear_modulus * scalar;
  const GenericReal<is_ad> creep_rate_derivative =
      -_coefficient * _three_shear_modulus * _n_exponent * pow(stress_delta, _n_exponent - 1.0) *
      _exponential * _exp_time;
  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
Real
PowerLawCreepStressUpdateTempl<is_ad>::computeStrainEnergyRateDensity(
    const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
    const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate)
{
  if (_n_exponent <= 1)
    return 0.0;

  Real creep_factor = _n_exponent / (_n_exponent + 1);

  return MetaPhysicL::raw_value(creep_factor * stress[_qp].doubleContraction((strain_rate)[_qp]));
}

template <bool is_ad>
void
PowerLawCreepStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
PowerLawCreepStressUpdateTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
PowerLawCreepStressUpdateTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->_use_substepping != RadialReturnStressUpdateTempl<is_ad>::SubsteppingType::NONE;
}

template class PowerLawCreepStressUpdateTempl<false>;
template class PowerLawCreepStressUpdateTempl<true>;
template Real PowerLawCreepStressUpdateTempl<false>::computeResidualInternal<Real>(const Real &,
                                                                                   const Real &);
template ADReal
PowerLawCreepStressUpdateTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                      const ADReal &);
template ChainedReal
PowerLawCreepStressUpdateTempl<false>::computeResidualInternal<ChainedReal>(const Real &,
                                                                            const ChainedReal &);
template ChainedADReal
PowerLawCreepStressUpdateTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &,
                                                                             const ChainedADReal &);
