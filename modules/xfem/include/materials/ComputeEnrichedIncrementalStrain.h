//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// #include "ComputeStrainBase.h"
#include "ComputeIncrementalStrainBase.h"
#include "ADComputeIncrementalStrainBase.h"
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "Assembly.h"
#include "CrackFrontDefinition.h"
#include "EnrichFunctionUtility.h"
#include "NonlinearSystem.h"

template <bool is_ad>
using ComputeIncrementalStrainBaseParent = typename std::
    conditional<is_ad, ADComputeIncrementalStrainBase, ComputeIncrementalStrainBase>::type;

template <bool is_ad>
using RealParent = typename std::conditional<is_ad, ADReal, Real>::type;

template <bool is_ad>
using RealVectorValueParent =
    typename std::conditional<is_ad, ADRealVectorValue, RealVectorValue>::type;

/**
 * ComputeIncrementalStrain defines a strain increment and rotation increment (=1), for small
 * strains.
 */
template <bool is_ad>
class ComputeEnrichedIncrementalStrainTempl : public ComputeIncrementalStrainBaseParent<is_ad>
{
public:
  static InputParameters validParams();

  ComputeEnrichedIncrementalStrainTempl(const InputParameters & parameters);
  virtual ~ComputeEnrichedIncrementalStrainTempl() {}
  virtual void initialSetup() override final;

  virtual void computeProperties() override;

protected:
  // /// Computes the current and old deformation gradients and passes back the
  // /// total strain increment tensor
  // virtual void computeTotalStrainIncrement(RankTwoTensor & total_strain_increment);

  /// enrichment displacement
  std::vector<RealParent<is_ad>> _enrich_disp;

  /// gradient of enrichment displacement
  std::vector<RealVectorValueParent<is_ad>> _grad_enrich_disp;
  std::vector<RealVectorValue> _grad_enrich_disp_old;

  /// enrichment displacement variables
  std::vector<std::vector<MooseVariableFEBase *>> _enrich_variable;

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// gradient of the shape function
  const VariablePhiGradient & _grad_phi;

  const CrackFrontDefinition * _crack_front_definition;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;
  GenericMaterialProperty<RankTwoTensor, is_ad> & _grad_disp_tensor;
  GenericMaterialProperty<RankTwoTensor, is_ad> & _small_strain;
  GenericMaterialProperty<RankTwoTensor, is_ad> & _grad_enrich_disp_tensor;
  const MaterialProperty<RankTwoTensor> & _grad_disp_tensor_old;
  const MaterialProperty<RankTwoTensor> & _small_strain_old;
  const MaterialProperty<RankTwoTensor> & _grad_enrich_disp_tensor_old;
  using ComputeIncrementalStrainBaseParent<is_ad>::_fe_problem;
  using ComputeIncrementalStrainBaseParent<is_ad>::_assembly;
  using ComputeIncrementalStrainBaseParent<is_ad>::_ndisp;
  using ComputeIncrementalStrainBaseParent<is_ad>::_current_elem;
  using ComputeIncrementalStrainBaseParent<is_ad>::_qrule;
  using ComputeIncrementalStrainBaseParent<is_ad>::isBoundaryMaterial;
  using ComputeIncrementalStrainBaseParent<is_ad>::_current_side;
  using ComputeIncrementalStrainBaseParent<is_ad>::_qp;
  using ComputeIncrementalStrainBaseParent<is_ad>::_q_point;
  using ComputeIncrementalStrainBaseParent<is_ad>::_grad_disp;
  // using ComputeIncrementalStrainBaseParent<is_ad>::_grad_phi;
  using ComputeIncrementalStrainBaseParent<is_ad>::_strain_increment;
  using ComputeIncrementalStrainBaseParent<is_ad>::_total_strain;
  using ComputeIncrementalStrainBaseParent<is_ad>::_dt;
  using ComputeIncrementalStrainBaseParent<is_ad>::_strain_rate;
  using ComputeIncrementalStrainBaseParent<is_ad>::_mechanical_strain;
  using ComputeIncrementalStrainBaseParent<is_ad>::_rotation_increment;

private:
  /// enrichment function value
  std::vector<Real> _B;
  /// derivatives of enrichment function respect to global cooridnate
  std::vector<RealVectorValue> _dBX;
  /// derivatives of enrichment function respect to crack front cooridnate
  std::vector<RealVectorValue> _dBx;
  /// enrichment function at node I
  std::vector<std::vector<Real>> _BI;
  /// shape function
  const std::vector<std::vector<Real>> * _fe_phi;
  /// gradient of shape function
  const std::vector<std::vector<RealGradient>> * _fe_dphi;
  NonlinearSystem * _nl;
  const NumericVector<Number> * _sln;

  using RealParentType = RealParent<is_ad>;
  using RealParentVectorType = RealVectorValueParent<is_ad>;
};
// typedef ComputeEnrichedIncrementalStrainTempl<true> ADComputeEnrichedIncrementalStrain;
// typedef ComputeEnrichedIncrementalStrainTempl<false> ComputeEnrichedIncrementalStrain;

using ComputeEnrichedIncrementalStrain = ComputeEnrichedIncrementalStrainTempl<false>;
using ADComputeEnrichedIncrementalStrain = ComputeEnrichedIncrementalStrainTempl<true>;
