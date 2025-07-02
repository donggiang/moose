//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeCrackTipEnrichmentIncrementalStrain.h"
#include "MooseMesh.h"
#include "libmesh/fe_interface.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/quadrature_gauss.h"

registerMooseObject("XFEMApp", ADComputeCrackTipEnrichmentIncrementalStrain);

InputParameters
ADComputeCrackTipEnrichmentIncrementalStrain::validParams()
{
  InputParameters params = ADComputeIncrementalStrainBase::validParams();
  params.addClassDescription(
      "Computes the crack tip enrichment at a point within a  incremental strain formulation.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>("enrichment_displacements",
                                                              "The enrichment displacement");
  params.addRequiredParam<UserObjectName>("crack_front_definition",
                                          "The CrackFrontDefinition user object name");
  return params;
}

ADComputeCrackTipEnrichmentIncrementalStrain::ADComputeCrackTipEnrichmentIncrementalStrain(
    const InputParameters & parameters)
  : ADComputeIncrementalStrainBase(parameters),
    _enrich_disp(3),
    _grad_enrich_disp(3),
    _grad_enrich_disp_old(3),
    _enrich_variable(4),
    _phi(_assembly.phi()),
    _grad_phi(_assembly.gradPhi()),
    _crack_front_definition(nullptr),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "total_strain")),
    _grad_disp_tensor(declareADProperty<RankTwoTensor>(_base_name + "grad_disp_tensor")),
    _small_strain(declareADProperty<RankTwoTensor>(_base_name + "small_strain")),
    _grad_enrich_disp_tensor(
        declareADProperty<RankTwoTensor>(_base_name + "grad_enrich_disp_tensor")),
    _grad_disp_tensor_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "grad_disp_tensor")),
    _small_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "small_strain")),
    _grad_enrich_disp_tensor_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "grad_enrich_disp_tensor")),
    _B(4),
    _dBX(4),
    _dBx(4)
{
  for (unsigned int i = 0; i < _enrich_variable.size(); ++i)
    _enrich_variable[i].resize(_ndisp);

  const std::vector<NonlinearVariableName> & nl_vnames =
      getParam<std::vector<NonlinearVariableName>>("enrichment_displacements");

  if (_ndisp == 2 && nl_vnames.size() != 8)
    mooseError("The number of enrichment displacements should be total 8 for 2D.");
  else if (_ndisp == 3 && nl_vnames.size() != 12)
    mooseError("The number of enrichment displacements should be total 12 for 3D.");

  _nl = &(_fe_problem.getNonlinearSystem(/*nl_sys_num=*/0));

  for (unsigned int j = 0; j < _ndisp; ++j)
    for (unsigned int i = 0; i < 4; ++i)
      _enrich_variable[i][j] = &(_nl->getVariable(0, nl_vnames[j * 4 + i]));

  if (_ndisp == 2)
    _BI.resize(4); // QUAD4
  else if (_ndisp == 3)
    _BI.resize(8); // HEX8

  for (unsigned int i = 0; i < _BI.size(); ++i)
    _BI[i].resize(4);
}

void
ADComputeCrackTipEnrichmentIncrementalStrain::initialSetup()
{
  const auto uo_name = getParam<UserObjectName>("crack_front_definition");
  _crack_front_definition = &_fe_problem.getUserObject<CrackFrontDefinition>(uo_name);
}

void
ADComputeCrackTipEnrichmentIncrementalStrain::computeProperties()
{
  // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug1\n";
  FEType fe_type(Utility::string_to_enum<Order>("first"),
                 Utility::string_to_enum<FEFamily>("lagrange"));
  const unsigned int dim = _current_elem->dim();
  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule(const_cast<QBase *>(_qrule));
  _fe_phi = &(fe->get_phi());
  _fe_dphi = &(fe->get_dphi());

  if (isBoundaryMaterial())
    fe->reinit(_current_elem, _current_side);
  else
    fe->reinit(_current_elem);

  _sln = _nl->currentSolution();

  for (unsigned int i = 0; i < _BI.size(); ++i)
    EnrichFunctionUtility::crackTipEnrichementFunctionAtPoint(
        _crack_front_definition, *(_current_elem->node_ptr(i)), _BI[i]);

  // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug2\n";

  // // fetch coupled gradients
  // for (unsigned int i = 0; i < _ndisp; ++i)
  //   _grad_disp[i] = &adCoupledGradients("displacements", i);

  // // fetch coupled gradients previous step
  // for (unsigned int i = 0; i < _ndisp; ++i)
  //   _grad_disp_old[i] = &coupledGradientOld("displacements", i);

  // incremental strain
  Real volumetric_strain = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    // enrichment function
    EnrichFunctionUtility::crackTipEnrichementFunctionAtPoint(
        _crack_front_definition, _q_point[_qp], _B);
    unsigned int crack_front_point_index =
        EnrichFunctionUtility::crackTipEnrichementFunctionDerivativeAtPoint(
            _crack_front_definition, _q_point[_qp], _dBx);
    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3\n";

    for (unsigned int i = 0; i < 4; ++i)
      EnrichFunctionUtility::rotateFromCrackFrontCoordsToGlobal(
          _crack_front_definition, _dBx[i], _dBX[i], crack_front_point_index);
    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.0\n";

    for (unsigned int m = 0; m < _ndisp; ++m)
    {
      _enrich_disp[m] = 0.0;
      _grad_enrich_disp[m].zero();
      // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.00\n";

      _grad_enrich_disp_old[m].zero();
      // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.000\n";

      for (unsigned int i = 0; i < _current_elem->n_nodes(); ++i)
      {
        // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.000\n";

        const Node * node_i = _current_elem->node_ptr(i);
        // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.0000\n";

        for (unsigned int j = 0; j < 4; ++j)
        {
          // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.1\n";

          dof_id_type dof = node_i->dof_number(_nl->number(), _enrich_variable[j][m]->number(), 0);
          ADReal soln = (*_sln)(dof);
          Real soln_old = _nl->solutionOld()(dof);
          // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.2\n";

          if (ADReal::do_derivatives)
            Moose::derivInsert(soln.derivatives(), dof, 1.);
          // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.3\n";

          _enrich_disp[m] += (*_fe_phi)[i][_qp] * (_B[j] - _BI[i][j]) * soln;
          RealVectorValue grad_B(_dBX[j]);
          // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.4\n";

          _grad_enrich_disp[m] +=
              ((*_fe_dphi)[i][_qp] * (_B[j] - _BI[i][j]) + (*_fe_phi)[i][_qp] * grad_B) * soln;
          // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug3.5\n";

          _grad_enrich_disp_old[m] +=
              ((*_fe_dphi)[i][_qp] * (_B[j] - _BI[i][j]) + (*_fe_phi)[i][_qp] * grad_B) * soln_old;
        }
      }
    }
    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug4\n";

    _grad_enrich_disp_tensor[_qp] = ADRankTwoTensor::initializeFromRows(
        _grad_enrich_disp[0], _grad_enrich_disp[1], _grad_enrich_disp[2]);
    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug4.1\n";

    RankTwoTensor grad_enrich_disp_tensor_old = RankTwoTensor::initializeFromRows(
        _grad_enrich_disp_old[0], _grad_enrich_disp_old[1], _grad_enrich_disp_old[2]);

    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug4.2\n";

    /////  standard strain part
    // Deformation gradient
    _grad_disp_tensor[_qp] = ADRankTwoTensor::initializeFromRows(
        (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug4.3\n";

    // auto grad_disp_tensor_old = RankTwoTensor ::initializeFromRows(
    //     (*_grad_disp_old[0])[_qp], (*_grad_disp_old[1])[_qp], (*_grad_disp_old[2])[_qp]);

    // _deformation_gradient[_qp] = _grad_disp_tensor[_qp] + _grad_enrich_disp_tensor[_qp];
    // _deformation_gradient[_qp].addIa(1.0);

    // std::cout << "ADComputeCrackTipEnrichmentIncrementalStrain, debug5\n";

    _small_strain[_qp] =
        0.5 * ((_grad_disp_tensor[_qp] + _grad_enrich_disp_tensor[_qp]) +
               (_grad_disp_tensor[_qp] + _grad_enrich_disp_tensor[_qp]).transpose());

    // RankTwoTensor small_strain_old =
    //     0.5 * ((grad_disp_tensor_old + grad_enrich_disp_tensor_old) +
    //            (grad_disp_tensor_old + grad_enrich_disp_tensor_old).transpose());

    _strain_increment[_qp] = _small_strain[_qp] - _small_strain_old[_qp];

    _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];

    // // Remove the Eigen strain increment
    // subtractEigenstrainIncrementFromStrain(_strain_increment[_qp]);

    // strain rate
    if (_dt > 0)
      _strain_rate[_qp] = _strain_increment[_qp] / _dt;
    else
      _strain_rate[_qp].zero();

    // Update strain in intermediate configuration: rotations are not needed
    _mechanical_strain[_qp] = _mechanical_strain_old[_qp] + _strain_increment[_qp];

    // incremental small strain does not include rotation
    _rotation_increment[_qp].setToIdentity();
  }
}
