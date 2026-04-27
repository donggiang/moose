//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef NEML2_ENABLED

#include <torch/script.h>
#include "neml2/models/Model.h"

namespace neml2
{

/**
 * NEML2 model that wraps a TorchScript Mixture-of-Experts surrogate (.pt) and
 * outputs only the **equivalent plastic strain rate** (ep_rate).
 *
 * Mirrors the role of LAROMANCE6DInterpolation with output_rate = ep_rate:
 * used both inside the J2 radial-return Newton solve and at the outer
 * composed-model level to provide the final ep_rate after convergence.
 *
 * For cell and wall dislocation density rates use MixtureOfExpertsCellWallModel.
 *
 * Inputs (physical units):
 *   von_mises_stress          [MPa]    effective stress
 *   temperature               [K]      absolute temperature
 *   equivalent_plastic_strain [-]      accumulated plastic strain
 *   cell_dislocation_density  [m^-2]   cell dislocation density
 *   wall_dislocation_density  [m^-2]   wall dislocation density
 *   flux                      [-]      normalised neutron flux in [0, 1]
 *
 * Output (physical units):
 *   ep_rate    [s^-1]        equivalent plastic strain rate
 */
class MixtureOfExpertsModel : public Model
{
public:
  static OptionSet expected_options();

  MixtureOfExpertsModel(const OptionSet & options);

  /// Send the TorchScript model to a different device / dtype.
  virtual void to(const torch::TensorOptions & options) override;

  void request_AD() override;

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

private:
  // ---- input variables (physical units) --------------------------------
  const Variable<Scalar> & _vm_stress;
  const Variable<Scalar> & _temperature;
  const Variable<Scalar> & _ep_strain;
  const Variable<Scalar> & _cell_dd;
  const Variable<Scalar> & _wall_dd;
  const Variable<Scalar> & _flux;

  // ---- output variable (ep_rate only) ----------------------------------
  Variable<Scalar> & _ep_rate;

  // ---- TorchScript surrogate -------------------------------------------
  std::unique_ptr<torch::jit::script::Module> _surrogate;
};

/**
 * NEML2 model that wraps a TorchScript Mixture-of-Experts surrogate (.pt) and
 * outputs **cell_rate and wall_rate** (dislocation density rates).
 *
 * Mirrors the combined role of LAROMANCE6DInterpolation with output_rate =
 * cell_rate / wall_rate: evaluated at the outer composed-model level after
 * the implicit Newton solve converges, using the converged state variables.
 *
 * Same inputs as MixtureOfExpertsModel; outputs two rates instead of one.
 */
class MixtureOfExpertsCellWallModel : public Model
{
public:
  static OptionSet expected_options();

  MixtureOfExpertsCellWallModel(const OptionSet & options);

  virtual void to(const torch::TensorOptions & options) override;

  void request_AD() override;

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

private:
  // ---- input variables (physical units) --------------------------------
  const Variable<Scalar> & _vm_stress;
  const Variable<Scalar> & _temperature;
  const Variable<Scalar> & _ep_strain;
  const Variable<Scalar> & _cell_dd;
  const Variable<Scalar> & _wall_dd;
  const Variable<Scalar> & _flux;

  // ---- output variables ------------------------------------------------
  Variable<Scalar> & _cell_rate;
  Variable<Scalar> & _wall_rate;

  // ---- TorchScript surrogate -------------------------------------------
  std::unique_ptr<torch::jit::script::Module> _surrogate;
};

} // namespace neml2

#endif // NEML2_ENABLED
