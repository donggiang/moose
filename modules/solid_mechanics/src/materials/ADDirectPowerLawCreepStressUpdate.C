//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ADDirectPowerLawCreepStressUpdate.h"

#include "RankTwoScalarTools.h"

registerMooseObject("SolidMechanicsApp", ADDirectPowerLawCreepStressUpdate);

InputParameters
ADDirectPowerLawCreepStressUpdate::validParams()
{
  InputParameters params = ADStressUpdateBase::validParams();
  params.addClassDescription(
      "Integrates isotropic power-law creep with a backward-Euler solve for the six stress "
      "components rather than a scalar radial-return solve.");
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addRequiredParam<Real>("coefficient", "Leading coefficient in the power-law equation");
  params.addRequiredRangeCheckedParam<Real>(
      "n_exponent", "n_exponent>=1", "Exponent on the von Mises stress");
  params.addParam<Real>("m_exponent", 0.0, "Exponent on time");
  params.addRequiredParam<Real>("activation_energy", "Activation energy");
  params.addParam<Real>("gas_constant", 8.3143, "Universal gas constant");
  params.addParam<Real>("start_time", 0.0, "Start time for the time-dependent factor");
  params.addRangeCheckedParam<Real>(
      "absolute_tolerance", 1e-12, "absolute_tolerance>0", "Local stress residual tolerance");
  params.addRangeCheckedParam<Real>(
      "relative_tolerance", 1e-10, "relative_tolerance>0", "Relative local tolerance");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations", 50, "max_iterations>0", "Maximum number of local Newton iterations");
  params.addParam<bool>("verbose", false, "Print the local Newton iteration count");
  params.addParamNamesToGroup("absolute_tolerance relative_tolerance max_iterations verbose",
                              "Advanced");
  return params;
}

ADDirectPowerLawCreepStressUpdate::ADDirectPowerLawCreepStressUpdate(
    const InputParameters & parameters)
  : ADStressUpdateBase(parameters),
    _temperature(isParamValid("temperature") ? &adCoupledValue("temperature") : nullptr),
    _coefficient(getParam<Real>("coefficient")),
    _n_exponent(getParam<Real>("n_exponent")),
    _m_exponent(getParam<Real>("m_exponent")),
    _activation_energy(getParam<Real>("activation_energy")),
    _gas_constant(getParam<Real>("gas_constant")),
    _start_time(getParam<Real>("start_time")),
    _absolute_tolerance(getParam<Real>("absolute_tolerance")),
    _relative_tolerance(getParam<Real>("relative_tolerance")),
    _max_iterations(getParam<unsigned int>("max_iterations")),
    _verbose(getParam<bool>("verbose")),
    _creep_strain(declareADProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _local_iterations(declareProperty<Real>(_base_name + "direct_stress_local_iterations"))
{
  if (_start_time < _app.getStartTime() && std::trunc(_m_exponent) != _m_exponent)
    paramError("start_time",
               "Start time must be equal to or greater than the Executioner start_time if a "
               "non-integer m_exponent is used");
}

void
ADDirectPowerLawCreepStressUpdate::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();
  _local_iterations[_qp] = 0.0;
}

void
ADDirectPowerLawCreepStressUpdate::propagateQpStatefulProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

ADRankTwoTensor
ADDirectPowerLawCreepStressUpdate::computeFlowRate(const ADRankTwoTensor & stress) const
{
  using std::exp;
  using std::pow;
  using std::sqrt;

  const ADRankTwoTensor deviatoric_stress = stress.deviatoric();
  const ADReal deviatoric_norm_squared = deviatoric_stress.doubleContraction(deviatoric_stress);
  if (MetaPhysicL::raw_value(deviatoric_norm_squared) == 0.0)
    return ADRankTwoTensor();

  const ADReal q = sqrt(1.5 * deviatoric_norm_squared);
  ADReal rate = _coefficient * pow(q, _n_exponent) * pow(_t - _start_time, _m_exponent);
  if (_temperature)
    rate *= exp(-_activation_energy / (_gas_constant * (*_temperature)[_qp]));

  return rate * 1.5 * deviatoric_stress / q;
}

void
ADDirectPowerLawCreepStressUpdate::computeResidualAndJacobian(
    const ADRankTwoTensor & stress,
    const ADRankTwoTensor & trial_stress,
    const ADRankFourTensor & elasticity_tensor,
    ADRankTwoTensor & residual,
    ADRankFourTensor & jacobian) const
{
  using std::pow;
  using std::sqrt;

  const ADRankTwoTensor deviatoric_stress = stress.deviatoric();
  const ADReal deviatoric_norm_squared = deviatoric_stress.doubleContraction(deviatoric_stress);

  residual = stress - trial_stress;
  jacobian = ADRankFourTensor(ADRankFourTensor::initIdentitySymmetricFour);

  if (MetaPhysicL::raw_value(deviatoric_norm_squared) == 0.0)
    return;

  const ADReal q = sqrt(1.5 * deviatoric_norm_squared);
  const ADRankTwoTensor flow_direction = 1.5 * deviatoric_stress / q;
  const ADRankTwoTensor flow_rate = computeFlowRate(stress);
  const ADReal factor = flow_rate.doubleContraction(flow_direction) /
                        flow_direction.doubleContraction(flow_direction);
  residual += _dt * (elasticity_tensor * flow_rate);

  const ADRankTwoTensor identity_two(ADRankTwoTensor::initIdentity);
  const ADRankFourTensor deviatoric_projection =
      ADRankFourTensor(ADRankFourTensor::initIdentitySymmetricFour) -
      identity_two.outerProduct(identity_two) / 3.0;
  const ADRankFourTensor flow_rate_derivative =
      (_n_exponent - 1.0) * factor / q * flow_direction.outerProduct(flow_direction) +
      1.5 * factor / q * deviatoric_projection;
  jacobian += _dt * elasticity_tensor * flow_rate_derivative;
}

void
ADDirectPowerLawCreepStressUpdate::computeVoigtJacobian(const ADRankFourTensor & jacobian,
                                                        DenseMatrix<ADReal> & voigt_jacobian) const
{
  for (const auto column : make_range(6))
  {
    ADRankTwoTensor basis;
    basis.zero();
    if (column < 3)
      basis(column, column) = 1.0;
    else
    {
      static const unsigned int first[] = {0, 1, 0};
      static const unsigned int second[] = {1, 2, 2};
      basis(first[column - 3], second[column - 3]) = 1.0;
      basis(second[column - 3], first[column - 3]) = 1.0;
    }

    const ADRankTwoTensor response = jacobian * basis;
    DenseVector<ADReal> response_vector(6);
    RankTwoScalarTools::RankTwoTensorToVoigtVector<true>(response, response_vector);
    for (const auto row : make_range(6))
      voigt_jacobian(row, column) = response_vector(row);
  }
}

void
ADDirectPowerLawCreepStressUpdate::updateState(ADRankTwoTensor & strain_increment,
                                               ADRankTwoTensor & inelastic_strain_increment,
                                               const ADRankTwoTensor & /*rotation_increment*/,
                                               ADRankTwoTensor & stress_new,
                                               const RankTwoTensor & /*stress_old*/,
                                               const ADRankFourTensor & elasticity_tensor,
                                               const RankTwoTensor & /*elastic_strain_old*/,
                                               bool /*compute_full_tangent_operator*/,
                                               RankFourTensor & /*tangent_operator*/)
{
  const ADRankTwoTensor trial_stress = stress_new;
  ADRankTwoTensor stress = trial_stress;
  unsigned int iteration = 0;
  bool converged = false;

  for (; iteration < _max_iterations; ++iteration)
  {
    ADRankTwoTensor residual;
    ADRankFourTensor jacobian;
    computeResidualAndJacobian(stress, trial_stress, elasticity_tensor, residual, jacobian);

    const Real residual_norm = MetaPhysicL::raw_value(residual.L2norm());
    const Real reference = std::max(MetaPhysicL::raw_value(trial_stress.L2norm()), 1.0);
    if (residual_norm <= _absolute_tolerance || residual_norm / reference <= _relative_tolerance)
    {
      converged = true;
      break;
    }

    DenseVector<ADReal> residual_vector(6);
    RankTwoScalarTools::RankTwoTensorToVoigtVector<true>(residual, residual_vector);
    residual_vector *= -1.0;

    DenseMatrix<ADReal> voigt_jacobian(6, 6);
    computeVoigtJacobian(jacobian, voigt_jacobian);
    DenseVector<ADReal> correction(6);
    voigt_jacobian.lu_solve(residual_vector, correction);

    ADRankTwoTensor stress_correction;
    RankTwoScalarTools::VoigtVectorToRankTwoTensor<true>(correction, stress_correction);
    stress += stress_correction;
  }

  if (!converged)
    mooseException("In ",
                   _name,
                   ": direct stress Newton solve did not converge in ",
                   _max_iterations,
                   " iterations.");

  _local_iterations[_qp] = iteration;
  if (_verbose)
    _console << _name << ": direct stress Newton converged in " << iteration << " iterations\n";

  inelastic_strain_increment = _dt * computeFlowRate(stress);
  _creep_strain[_qp] = _creep_strain_old[_qp] + inelastic_strain_increment;
  strain_increment -= inelastic_strain_increment;
  stress_new = stress;
}
