//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MixtureOfExpertsModel.h"

#ifdef NEML2_ENABLED

#include <torch/script.h>
#include <torch/torch.h>
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/abs.h"
#include "neml2/tensors/functions/sign.h"

namespace neml2
{
register_NEML2_object(MixtureOfExpertsModel);
register_NEML2_object(MixtureOfExpertsCellWallModel);

// ---------------------------------------------------------------------------
// Shared helper: declare the six standard input options
// ---------------------------------------------------------------------------
static void
add_moe_input_options(OptionSet & options)
{
  options.set<std::string>("file_path");
  options.set("file_path").doc() = "Path to the TorchScript (.pt) model file";

  options.set<bool>("jit") = false;
  options.set("jit").suppressed() = true;

  options.set<VariableName>("von_mises_stress");
  options.set("von_mises_stress").doc() = "von Mises effective stress [MPa]";

  options.set<VariableName>("temperature");
  options.set("temperature").doc() = "Absolute temperature [K]";

  options.set<VariableName>("equivalent_plastic_strain");
  options.set("equivalent_plastic_strain").doc() = "Accumulated equivalent plastic strain [-]";

  options.set<VariableName>("cell_dislocation_density");
  options.set("cell_dislocation_density").doc() = "Cell dislocation density [m^-2]";

  options.set<VariableName>("wall_dislocation_density");
  options.set("wall_dislocation_density").doc() = "Wall dislocation density [m^-2]";

  options.set<VariableName>("flux");
  options.set("flux").doc() = "Normalised neutron flux in [0, 1]";
}

// ---------------------------------------------------------------------------
// Shared helper: build [N, 6] normalised input tensor from the six inputs
// ---------------------------------------------------------------------------
static Tensor
build_moe_input(const Variable<Scalar> & vm_stress,
                const Variable<Scalar> & temperature,
                const Variable<Scalar> & ep_strain,
                const Variable<Scalar> & cell_dd,
                const Variable<Scalar> & wall_dd,
                const Variable<Scalar> & flux)
{
  const auto n_batch_dims = vm_stress.batch_dim();

  // col 0: von Mises stress [MPa], minmax
  constexpr double s_min = 0.11428862639259419753390289997697;
  constexpr double s_max = 299.94931150862146296276478096842766;
  const auto v1 = (vm_stress() - s_min) / ((s_max - s_min) / 2.0) - 1.0;

  // col 1: temperature [K], minmax
  constexpr double t_min = 600.04429302445396388066001236438751;
  constexpr double t_max = 1099.98408980407225499220658093690872;
  const auto v2 = (temperature() - t_min) / ((t_max - t_min) / 2.0) - 1.0;

  // col 2: equivalent plastic strain [-], log-transform (eps = 1e-20)
  constexpr double ep_logmin = -39.15241605429763183110480895265937;
  constexpr double ep_logmax = -3.91206356603903770974284270778298;
  const auto v3 =
      (neml2::log(ep_strain() + 1e-20) - ep_logmin) / ((ep_logmax - ep_logmin) / 2.0) - 1.0;

  // col 3: cell dislocation density [m^-2], log-transform
  constexpr double c_logmin = 12.56521587815730534032354626106098;
  constexpr double c_logmax = 29.64417981061014728538793860934675;
  const auto v4 = (neml2::log(cell_dd()) - c_logmin) / ((c_logmax - c_logmin) / 2.0) - 1.0;

  // col 4: wall dislocation density [m^-2], log-transform
  constexpr double w_logmin = 29.12196678204971789227784029208124;
  constexpr double w_logmax = 30.11574121515777946456182689871639;
  const auto v5 = (neml2::log(wall_dd()) - w_logmin) / ((w_logmax - w_logmin) / 2.0) - 1.0;

  // col 5: neutron flux [-], minmax (xmin=0, xmax=1)
  const auto v6 = flux() * 2.0 - 1.0;

  const std::vector<at::Tensor> cols = {static_cast<at::Tensor>(v1),
                                        static_cast<at::Tensor>(v2),
                                        static_cast<at::Tensor>(v3),
                                        static_cast<at::Tensor>(v4),
                                        static_cast<at::Tensor>(v5),
                                        static_cast<at::Tensor>(v6)};
  return Tensor(torch::transpose(torch::vstack(cols), 0, 1), n_batch_dims);
}

// ===========================================================================
// MixtureOfExpertsModel  (ep_rate only)
// ===========================================================================

OptionSet
MixtureOfExpertsModel::expected_options()
{
  auto options = Model::expected_options();
  options.doc() =
      "Wraps a TorchScript Mixture-of-Experts surrogate (.pt) and outputs the "
      "equivalent plastic strain rate (ep_rate).  Apply input normalization and "
      "output de-normalization internally so the raw .pt file can be used directly.";

  add_moe_input_options(options);

  options.set_output("ep_rate");
  options.set("ep_rate").doc() = "Equivalent plastic strain rate [s^-1]";

  return options;
}

MixtureOfExpertsModel::MixtureOfExpertsModel(const OptionSet & options)
  : Model(options),
    _vm_stress(declare_input_variable<Scalar>("von_mises_stress")),
    _temperature(declare_input_variable<Scalar>("temperature")),
    _ep_strain(declare_input_variable<Scalar>("equivalent_plastic_strain")),
    _cell_dd(declare_input_variable<Scalar>("cell_dislocation_density")),
    _wall_dd(declare_input_variable<Scalar>("wall_dislocation_density")),
    _flux(declare_input_variable<Scalar>("flux")),
    _ep_rate(declare_output_variable<Scalar>("ep_rate")),
    _surrogate(std::make_unique<torch::jit::script::Module>(
        torch::jit::load(options.get<std::string>("file_path"))))
{
  // Force single-threaded LibTorch evaluation so results are independent of
  // how many MPI processes split the Gauss-point batch (parallel reductions
  // are non-associative in floating-point; different thread counts → different
  // rounding order → irreproducible results across ranks).
  // Force single intra-op thread so floating-point reductions inside LibTorch
  // tensor ops are evaluated in a fixed order, making results independent of
  // how many MPI processes split the Gauss-point batch.
  at::set_num_threads(1);
}

void
MixtureOfExpertsModel::to(const torch::TensorOptions & options)
{
  Model::to(options);
  if (options.has_device())
    _surrogate->to(options.device());
  if (options.has_dtype())
    _surrogate->to(torch::Dtype(caffe2::typeMetaToScalarType(options.dtype())));
}

void
MixtureOfExpertsModel::request_AD()
{
  const std::vector<const VariableBase *> inputs = {
      &_vm_stress, &_temperature, &_ep_strain, &_cell_dd, &_wall_dd, &_flux};
  _ep_rate.request_AD(inputs);
}

void
MixtureOfExpertsModel::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto x_norm =
      build_moe_input(_vm_stress, _temperature, _ep_strain, _cell_dd, _wall_dd, _flux);

  const auto raw = _surrogate->forward({x_norm}).toTensor(); // [N, 3]

  // out col 0: ep_rate [s^-1], logtransform_inv (lb=0, ub=1)
  constexpr double lxmin = -30.54946341410551013950680498965085;
  constexpr double lxmax = 7.72739242062932873977842973545194;
  const Scalar z = Scalar(raw.select(1, 0), 0);
  _ep_rate = neml2::exp(z * (lxmax - lxmin) + lxmin);
}

// ===========================================================================
// MixtureOfExpertsCellWallModel  (cell_rate and wall_rate)
// ===========================================================================

OptionSet
MixtureOfExpertsCellWallModel::expected_options()
{
  auto options = Model::expected_options();
  options.doc() =
      "Wraps a TorchScript Mixture-of-Experts surrogate (.pt) and outputs the "
      "cell and wall dislocation density rates (cell_rate, wall_rate).  Apply "
      "input normalization and output de-normalization internally so the raw .pt "
      "file can be used directly.";

  add_moe_input_options(options);

  options.set_output("cell_rate");
  options.set("cell_rate").doc() = "Cell dislocation density rate [m^-2 s^-1]";

  options.set_output("wall_rate");
  options.set("wall_rate").doc() = "Wall dislocation density rate [m^-2 s^-1]";

  return options;
}

MixtureOfExpertsCellWallModel::MixtureOfExpertsCellWallModel(const OptionSet & options)
  : Model(options),
    _vm_stress(declare_input_variable<Scalar>("von_mises_stress")),
    _temperature(declare_input_variable<Scalar>("temperature")),
    _ep_strain(declare_input_variable<Scalar>("equivalent_plastic_strain")),
    _cell_dd(declare_input_variable<Scalar>("cell_dislocation_density")),
    _wall_dd(declare_input_variable<Scalar>("wall_dislocation_density")),
    _flux(declare_input_variable<Scalar>("flux")),
    _cell_rate(declare_output_variable<Scalar>("cell_rate")),
    _wall_rate(declare_output_variable<Scalar>("wall_rate")),
    _surrogate(std::make_unique<torch::jit::script::Module>(
        torch::jit::load(options.get<std::string>("file_path"))))
{
  // Force single intra-op thread so floating-point reductions inside LibTorch
  // tensor ops are evaluated in a fixed order, making results independent of
  // how many MPI processes split the Gauss-point batch.
  at::set_num_threads(1);
}

void
MixtureOfExpertsCellWallModel::to(const torch::TensorOptions & options)
{
  Model::to(options);
  if (options.has_device())
    _surrogate->to(options.device());
  if (options.has_dtype())
    _surrogate->to(torch::Dtype(caffe2::typeMetaToScalarType(options.dtype())));
}

void
MixtureOfExpertsCellWallModel::request_AD()
{
  const std::vector<const VariableBase *> inputs = {
      &_vm_stress, &_temperature, &_ep_strain, &_cell_dd, &_wall_dd, &_flux};
  _cell_rate.request_AD(inputs);
  _wall_rate.request_AD(inputs);
}

void
MixtureOfExpertsCellWallModel::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto x_norm =
      build_moe_input(_vm_stress, _temperature, _ep_strain, _cell_dd, _wall_dd, _flux);

  const auto raw = _surrogate->forward({x_norm}).toTensor(); // [N, 3]

  // out col 1: cell_rate [m^-2 s^-1], symlog_inv
  {
    constexpr double zmin = -36.71520922998040958873389172367752;
    constexpr double zmax = 38.45450753607088500984900747425854;
    const Scalar zs = Scalar(raw.select(1, 1), 0);
    const auto z = (zs + 1.0) * ((zmax - zmin) / 2.0) + zmin;
    _cell_rate = neml2::sign(z) * (neml2::exp(neml2::abs(z)) - 1.0);
  }

  // out col 2: wall_rate [m^-2 s^-1], symlog_inv
  {
    constexpr double zmin = -40.36444506932617315442257677204907;
    constexpr double zmax = 35.04333232507948281408971524797380;
    const Scalar zs = Scalar(raw.select(1, 2), 0);
    const auto z = (zs + 1.0) * ((zmax - zmin) / 2.0) + zmin;
    _wall_rate = neml2::sign(z) * (neml2::exp(neml2::abs(z)) - 1.0);
  }
}

} // namespace neml2

#endif // NEML2_ENABLED
