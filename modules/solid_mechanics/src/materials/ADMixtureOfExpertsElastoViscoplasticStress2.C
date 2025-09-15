//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMixtureOfExpertsElastoViscoplasticStress2.h"
#include <torch/script.h> // One-stop header.
#include <torch/torch.h>
#include <iostream>
#include <memory>

registerMooseObject("SolidMechanicsApp", ADMixtureOfExpertsElastoViscoplasticStress2);

InputParameters
ADMixtureOfExpertsElastoViscoplasticStress2::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription(
      "Computes stress after subtracting inelastic strain increment calculated with"
      "Mixture of Expert ML based elasto-viscoplastic constitutive relation");
  params.addRequiredParam<MaterialPropertyName>("mu", "Shear modulus");
  params.addRequiredParam<MaterialPropertyName>("temp", "Temperature");
  params.addRequiredParam<Real>("initial_evm", "Initial homogeneous effective strain");
  params.addRequiredParam<Real>("initial_rhoc", "Initial Dislocation density at cell");
  params.addRequiredParam<Real>("initial_rhow", "Initial Dislocation density at cell wall");
  params.addRequiredParam<MaterialPropertyName>("flux", "Neutron Flux");
  return params;
}

ADMixtureOfExpertsElastoViscoplasticStress2::ADMixtureOfExpertsElastoViscoplasticStress2(
    const InputParameters & parameters)
  : ADComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _mechanical_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mu(getADMaterialProperty<Real>("mu")),
    _temp(getADMaterialProperty<Real>("temp")),
    _rhoc(declareADProperty<Real>("rhoc")),
    _rhow(declareADProperty<Real>("rhow")),
    _rhoc_old(getMaterialPropertyOld<Real>("rhoc")),
    _rhow_old(getMaterialPropertyOld<Real>("rhow")),
    _flux(getADMaterialProperty<Real>("flux")),
    _total_effective_plastic_strain(declareADProperty<Real>("effective_plastic_strain")),
    _total_effective_plastic_strain_old(getMaterialPropertyOld<Real>("effective_plastic_strain")),
    _effective_plastic_strain_rate(declareADProperty<Real>("effective_plastic_strain_rate")),
    _initial_evm(getParam<Real>("initial_evm")),
    _initial_rhoc(getParam<Real>("initial_rhoc")),
    _initial_rhow(getParam<Real>("initial_rhow")),
    effective_plastic_strain_increment(0.0),
    rhoc_rate(0.0),
    rhow_rate(0.0),
    _incremental_plastic_strain_tensor(
        declareADProperty<RankTwoTensor>("incremental_plastic_strain_tensor")),
    _incremental_plastic_strain_tensor_old(
        getMaterialPropertyOld<RankTwoTensor>("incremental_plastic_strain_tensor"))

{
  // Checking if the model file exists and load it during construction
  std::string model_path = "/Users/huyngd/projects/moose-torch/model9experts_new.pt";
  // std::string model_path = "/Users/zaheennasir/projects/snail/model2experts.pt";

  if (!std::filesystem::exists(model_path))
  {
    std::cerr << "Model file does not exist at the specified path: " << model_path << std::endl;
    throw std::runtime_error("Model file not found");
  }

  try
  {
    // Load the TorchScript model
    std::cout << "Stess: " << std::endl;
    model = torch::jit::load(model_path);
    std::cout << "Model loaded successfully." << std::endl;
  }
  catch (const c10::Error & e)
  {
    std::cerr << "Error loading the model: " << e.what() << std::endl;
    throw; // Propagate the exception
  }
}

void
ADMixtureOfExpertsElastoViscoplasticStress2::initQpStatefulProperties()
{
  _total_effective_plastic_strain[_qp] = _initial_evm;
  _rhoc[_qp] = _initial_rhoc;
  _rhow[_qp] = _initial_rhow;
}

void
ADMixtureOfExpertsElastoViscoplasticStress2::computeQpStress()
{
  // copying dislocations to a temporary file so that any change can be enforced
  auto _rhoc_old_temp = _rhoc_old[_qp];
  auto _rhow_old_temp = _rhow_old[_qp];

  // Checking bounds of the MoE model inputs and fixing them
  if (_rhoc_old_temp <= 0)
  {
    _rhoc_old_temp = _initial_rhoc;
  }

  if (_rhow_old_temp <= 0)
  {
    _rhow_old_temp = _initial_rhow;
  }

  if (_rhoc_old_temp > 8461801410123.313)
  {
    _rhoc_old_temp = 8461801410123.313;
  }

  if (_rhow_old_temp > 11999567054170.322)
  {
    _rhow_old_temp = 11999567054170.322;
  }

  // Getting the strain increment
  ADRankTwoTensor elastic_strain_increment =
      (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]);

  // Computing trial stress asuming the whole strain increment as elastic
  RankTwoTensor trial_stress =
      MetaPhysicL::raw_value(_stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment);
  ADRankTwoTensor deviatoric_trial_stress = trial_stress.deviatoric();
  auto norm_dev_stress_squared = deviatoric_trial_stress.doubleContraction(deviatoric_trial_stress);
  auto norm_dev_stress = std::sqrt(norm_dev_stress_squared);
  auto effective_trial_stress = std::sqrt(1.5) * norm_dev_stress;
  const Real tolerance = 1e-14;

  // Checking for non-zero stress

  if (effective_trial_stress > (1.0 + tolerance))
  {

    // Creating an input tensor
    at::Tensor inputs =
        torch::tensor({{MetaPhysicL::raw_value(effective_trial_stress),
                        MetaPhysicL::raw_value(_temp[_qp]),
                        MetaPhysicL::raw_value(_total_effective_plastic_strain_old[_qp]),
                        MetaPhysicL::raw_value(_rhoc_old_temp),
                        MetaPhysicL::raw_value(_rhow_old_temp),
                        MetaPhysicL::raw_value(_flux[_qp])}},
                      torch::kFloat);

    // evaluating the model
    auto model_output = this->eval(inputs);
    auto creep_rate = model_output.index({0, 0});
    // initializing useful varialbles
    Real effective_plastic_strain_increment = 0.0;
    Real residual = 0.0;
    unsigned int k = 0;
    Real tol = 1e-8;
    Real dt = _dt;
    if (dt == 0)
    {
      dt = 1;
    }

    // Starting the Return-Mapping Newton-Raphson loop
    while (std::abs(residual) > tol || k == 0)
    {
      if (++k > 1000)
      {
        std::cerr << "Max iterations reached.\n";
        break;
      }

      residual = effective_plastic_strain_increment - creep_rate.item().toDouble() * dt;
      auto jacobian = this->eval_jacobian(inputs);
      auto jac_1st = jacobian.index({0, 0});
      Real dphi_ddeqpl = -3.0 * MetaPhysicL::raw_value(_mu[_qp]) * jac_1st.item().toDouble();
      Real del_deqpl = (creep_rate.item().toDouble() - effective_plastic_strain_increment / dt) /
                       ((1 / dt) - dphi_ddeqpl);
      effective_plastic_strain_increment += del_deqpl;
      at::Tensor inputs_new = torch::tensor(
          {{MetaPhysicL::raw_value(effective_trial_stress -
                                   3.0 * _mu[_qp] * effective_plastic_strain_increment),
            MetaPhysicL::raw_value(_temp[_qp]),
            MetaPhysicL::raw_value(_total_effective_plastic_strain_old[_qp]),
            MetaPhysicL::raw_value(_rhoc_old_temp),
            MetaPhysicL::raw_value(_rhow_old_temp),
            MetaPhysicL::raw_value(_flux[_qp])}},
          torch::kFloat);
      inputs = inputs_new;
      model_output = this->eval(inputs);
      auto new_rate = model_output.index({0, 0});
      creep_rate = new_rate;
    }

    // Updating dislocations and total plastic strains for next MoE input
    auto rhoc_rate_temp = model_output.index({0, 1});
    auto rhow_rate_temp = model_output.index({0, 2});

    rhoc_rate = rhoc_rate_temp.item().toDouble();
    rhow_rate = rhow_rate_temp.item().toDouble();

    _total_effective_plastic_strain[_qp] =
        _total_effective_plastic_strain_old[_qp] + effective_plastic_strain_increment;
    _rhoc[_qp] = _rhoc_old_temp + rhoc_rate * dt;
    _rhow[_qp] = _rhow_old_temp + rhow_rate * dt;
    _effective_plastic_strain_rate[_qp] = effective_plastic_strain_increment / dt;

    // Getting the plastic strain tensor
    _incremental_plastic_strain_tensor[_qp] =
        deviatoric_trial_stress *
        (1.5 * effective_plastic_strain_increment / effective_trial_stress);

    // getting elastic strain
    elastic_strain_increment = elastic_strain_increment - _incremental_plastic_strain_tensor[_qp];
  }
  // Computing actual stress after returning back
  _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::eval(const torch::Tensor & tensor)
{

  auto varx1 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 0}),
      -1.0,
      1.0,
      0.11428862639259419753390289997697,
      299.94931150862146296276478096842766);
  auto varx2 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 1}),
      -1.0,
      1.0,
      600.04429302445396388066001236438751,
      1099.98408980407225499220658093690872);
  auto varx3 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 2}),
      -1.0,
      1.0,
      -39.15241605429763183110480895265937,
      -3.91206356603903770974284270778298,
      1e-20);
  auto varx4 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 3}),
      -1.0,
      1.0,
      12.56521587815730534032354626106098,
      29.64417981061014728538793860934675,
      0.0);
  auto varx5 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 4}),
      -1.0,
      1.0,
      29.12196678204971789227784029208124,
      30.11574121515777946456182689871639,
      0.0);
  auto varx6 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  at::Tensor transformed = torch::tensor({{varx1.item<double>(),
                                           varx2.item<double>(),
                                           varx3.item<double>(),
                                           varx4.item<double>(),
                                           varx5.item<double>(),
                                           varx6.item<double>()}},
                                         torch::kFloat);

  at::Tensor output = model.forward({transformed}).toTensor();

  auto out1 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_inverse_transform(
      output.index({0, 0}),
      0.0,
      1.0,
      -30.54946341410551013950680498965085,
      7.72739242062932873977842973545194,
      0.0);
  auto out2 = ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_transform(
      output.index({0, 1}),
      -36.71520922998040958873389172367752,
      38.45450753607088500984900747425854);
  auto out3 = ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_transform(
      output.index({0, 2}),
      -40.36444506932617315442257677204907,
      35.04333232507948281408971524797380);

  at::Tensor evaluated = torch::tensor(
      {{out1.item<double>(), out2.item<double>(), out3.item<double>()}}, torch::kFloat);

  return evaluated;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::eval_jacobian(const torch::Tensor & tensor)
{
  auto jax1 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_forward_derivative(
      tensor.index({0, 0}),
      -1.0,
      1.0,
      0.11428862639259419753390289997697,
      299.94931150862146296276478096842766);
  auto jax2 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_forward_derivative(
      tensor.index({0, 1}),
      -1.0,
      1.0,
      600.04429302445396388066001236438751,
      1099.98408980407225499220658093690872);
  auto jax3 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_forward_derivative(
      tensor.index({0, 2}),
      -1.0,
      1.0,
      -39.15241605429763183110480895265937,
      -3.91206356603903770974284270778298,
      1e-20);
  auto jax4 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_forward_derivative(
      tensor.index({0, 3}),
      -1.0,
      1.0,
      12.56521587815730534032354626106098,
      29.64417981061014728538793860934675,
      0.0);
  auto jax5 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_forward_derivative(
      tensor.index({0, 4}),
      -1.0,
      1.0,
      29.12196678204971789227784029208124,
      30.11574121515777946456182689871639,
      0.0);
  auto jax6 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_forward_derivative(
      tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  at::Tensor Jin = torch::tensor({{jax1.item<double>(),
                                   jax2.item<double>(),
                                   jax3.item<double>(),
                                   jax4.item<double>(),
                                   jax5.item<double>(),
                                   jax6.item<double>()}},
                                 torch::kFloat)
                       .view(-1);

  auto varx1 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 0}),
      -1.0,
      1.0,
      0.11428862639259419753390289997697,
      299.94931150862146296276478096842766);
  auto varx2 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 1}),
      -1.0,
      1.0,
      600.04429302445396388066001236438751,
      1099.98408980407225499220658093690872);
  auto varx3 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 2}),
      -1.0,
      1.0,
      -39.15241605429763183110480895265937,
      -3.91206356603903770974284270778298,
      1e-20);
  auto varx4 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 3}),
      -1.0,
      1.0,
      12.56521587815730534032354626106098,
      29.64417981061014728538793860934675,
      0.0);
  auto varx5 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
      tensor.index({0, 4}),
      -1.0,
      1.0,
      29.12196678204971789227784029208124,
      30.11574121515777946456182689871639,
      0.0);
  auto varx6 = ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
      tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  at::Tensor p = torch::tensor({{varx1.item<double>(),
                                 varx2.item<double>(),
                                 varx3.item<double>(),
                                 varx4.item<double>(),
                                 varx5.item<double>(),
                                 varx6.item<double>()}},
                               torch::kFloat);

  p = p.repeat({3, 1});
  p.requires_grad_(true);

  torch::Tensor y = model.forward({p}).toTensor();

  torch::Tensor eye_matrix = torch::eye(3, torch::kDouble);

  y.backward(eye_matrix);
  torch::Tensor Jmod = p.grad();

  auto phi = y[0].unsqueeze(0);

  auto out1 = ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_inverse_derivative(
      phi.index({0, 0}),
      0.0,
      1.0,
      -30.54946341410551013950680498965085,
      7.72739242062932873977842973545194,
      0.0);
  auto out2 = ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_derivative(
      phi.index({0, 1}), -36.71520922998040958873389172367752, 38.45450753607088500984900747425854);
  auto out3 = ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_derivative(
      phi.index({0, 2}), -40.36444506932617315442257677204907, 35.04333232507948281408971524797380);

  at::Tensor Jout = torch::tensor({{out1.item<double>(), out2.item<double>(), out3.item<double>()}},
                                  torch::kFloat)
                        .view(-1);

  torch::Tensor J = Jout.unsqueeze(1) * (Jmod * Jin);

  return J;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::minmax_transform(
    const torch::Tensor & x, double lb, double ub, double xmin, double xmax)
{

  auto a = ub - lb;
  auto b = lb;
  return a * (x - xmin) / (xmax - xmin) + b;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_transform(
    const torch::Tensor & x, double lb, double ub, double logxmin, double logxmax, double eps)
{
  auto a = ub - lb;
  auto b = lb;
  return a * (torch::log(x + eps) - logxmin) / (logxmax - logxmin) + b;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_inverse_transform(
    const torch::Tensor & z, double lb, double ub, double logxmin, double logxmax, double eps)
{
  auto a = ub - lb;
  auto b = lb;
  return torch::exp((logxmax - logxmin) * (z - b) / a + logxmin) - eps;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_transform(
    const torch::Tensor & z_scaled, double zmin, double zmax)
{
  auto z = (z_scaled + 1) * (zmax - zmin) / 2 + zmin;

  // Then inverse symlog
  return torch::sign(z) * (torch::exp(torch::abs(z)) - 1.0);
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::minmax_forward_derivative(
    const torch::Tensor & x, double lb, double ub, double xmin, double xmax)
{
  auto a = ub - lb;
  // auto b = lb;
  return torch::tensor(a / (xmax - xmin));
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_forward_derivative(
    const torch::Tensor & x, double lb, double ub, double logxmin, double logxmax, double eps)
{
  auto a = ub - lb;
  auto b = lb;
  return a / (logxmax - logxmin) * (1.0 / ((x + eps)));
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::logtransform_inverse_derivative(
    const torch::Tensor & z, double lb, double ub, double logxmin, double logxmax, double eps)
{
  auto a = ub - lb;
  auto b = lb;
  auto u = torch::exp((logxmax - logxmin) * (z - b) / a + logxmin);
  return u * (logxmax - logxmin) / a;
}

torch::Tensor
ADMixtureOfExpertsElastoViscoplasticStress2::symlog_inverse_derivative(
    const torch::Tensor & z_scaled, double zmin, double zmax)
{
  auto z = (z_scaled + 1) * (zmax - zmin) / 2 + zmin;
  auto u = torch::exp(torch::abs(z));
  return u * (zmax - zmin) / 2;
}
