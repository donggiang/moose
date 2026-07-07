//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshCut2DCCGUserObject.h"

#include "CrackFrontDefinition.h"
#include "MooseError.h"

#include <cmath>

registerMooseObject("XFEMApp", MeshCut2DCCGUserObject);

InputParameters
MeshCut2DCCGUserObject::validParams()
{
  InputParameters params = MeshCut2DUserObjectBase::validParams();
  params.addClassDescription("XFEM mesh cutter for 2D models that defines cuts with a mesh and "
                             "uses the C integral to determine creep crack growth.");
  params.addRequiredParam<Real>("growth_increment",
                                "Retained for input compatibility. CCG growth is computed from "
                                "the C integral and Paris-law coefficients.");
  params.addRequiredParam<Real>("paris_coeff", "Coefficient in Paris's law.");
  params.addRequiredParam<Real>("paris_exponent", "Exponent in Paris's law.");
  params.addParam<VectorPostprocessorName>(
      "ki_vectorpostprocessor", "II_KI_1", "The name of the vectorpostprocessor that contains KI");
  params.addParam<VectorPostprocessorName>("kii_vectorpostprocessor",
                                           "II_KII_1",
                                           "The name of the vectorpostprocessor that contains KII");
  params.addParam<VectorPostprocessorName>(
      "c_vectorpostprocessor", "C_1", "The name of the vectorpostprocessor that contains C");
  return params;
}

MeshCut2DCCGUserObject::MeshCut2DCCGUserObject(const InputParameters & parameters)
  : MeshCut2DUserObjectBase(parameters),
    _ki_vpp(getVectorPostprocessorValue(
        "ki_vectorpostprocessor", getParam<VectorPostprocessorName>("ki_vectorpostprocessor"))),
    _kii_vpp(getVectorPostprocessorValue(
        "kii_vectorpostprocessor", getParam<VectorPostprocessorName>("kii_vectorpostprocessor"))),
    _c_vpp(getVectorPostprocessorValue("c_vectorpostprocessor",
                                       getParam<VectorPostprocessorName>("c_vectorpostprocessor"))),
    _total_crack_length(0.0),
    _paris_coeff(getParam<Real>("paris_coeff")),
    _paris_exponent(getParam<Real>("paris_exponent"))
{
}

void
MeshCut2DCCGUserObject::initialize()
{
  _is_mesh_modified = false;
  findActiveBoundaryGrowth();
  growFront();
  addNucleatedCracksToMesh();
  _crack_front_definition->updateNumberOfCrackFrontPoints(
      _original_and_current_front_node_ids.size());
  _crack_front_definition->isCutterModified(_is_mesh_modified);
  if (_is_mesh_modified)
    _crack_front_definition->updateCrackFrontPoints();
}

bool
MeshCut2DCCGUserObject::isCutterMeshChanged() const
{
  return _is_mesh_modified;
}

void
MeshCut2DCCGUserObject::findActiveBoundaryGrowth()
{
  // The C-integral vector is empty on the first XFEM update because this user object can execute
  // before the DomainIntegral vectorpostprocessor has produced values.
  if (_c_vpp.empty())
    return;

  if (_ki_vpp.size() != _kii_vpp.size() ||
      _ki_vpp.size() != _original_and_current_front_node_ids.size() ||
      _c_vpp.size() != _original_and_current_front_node_ids.size())
    mooseError("ki_vectorpostprocessor, kii_vectorpostprocessor, and c_vectorpostprocessor should "
               "have the same number of crack tips as CrackFrontDefinition.",
               "\n  ki size = ",
               _ki_vpp.size(),
               "\n  kii size = ",
               _kii_vpp.size(),
               "\n  C size = ",
               _c_vpp.size(),
               "\n  cracktips in MeshCut2DCCGUserObject = ",
               _original_and_current_front_node_ids.size());

  _active_front_node_growth_vectors.clear();
  for (const auto i : make_range(_original_and_current_front_node_ids.size()))
  {
    const Real c = _c_vpp[i];
    if (c <= 0.0)
      continue;

    const Real growth_increment = _paris_coeff * std::pow(c, _paris_exponent) * _dt;
    if (growth_increment <= 0.0)
      continue;

    _total_crack_length += growth_increment;

    const Real ki = _ki_vpp[i];
    const Real kii = _kii_vpp[i];
    const Real sqrt_k = std::sqrt(ki * ki + 8.0 * kii * kii);

    Real theta_m = 0.0;
    Real theta_p = 0.0;
    if (std::abs(kii) > libMesh::TOLERANCE)
    {
      theta_m = 2.0 * std::atan((ki - sqrt_k) / (4.0 * kii));
      theta_p = 2.0 * std::atan((ki + sqrt_k) / (4.0 * kii));
    }

    const Real sigma_tt_m = ki * (3.0 * std::cos(theta_m / 2.0) + std::cos(3.0 * theta_m / 2.0)) +
                            kii * (-3.0 * std::sin(theta_m / 2.0) -
                                   3.0 * std::sin(3.0 * theta_m / 2.0));
    const Real sigma_tt_p = ki * (3.0 * std::cos(theta_p / 2.0) + std::cos(3.0 * theta_p / 2.0)) +
                            kii * (-3.0 * std::sin(theta_p / 2.0) -
                                   3.0 * std::sin(3.0 * theta_p / 2.0));
    const Real theta = sigma_tt_m > sigma_tt_p ? theta_m : theta_p;

    RealVectorValue dir_cfc;
    dir_cfc(0) = std::cos(theta);
    dir_cfc(1) = std::sin(theta);
    dir_cfc(2) = 0.0;

    const RealVectorValue dir_global =
        _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);
    const Point nodal_offset(dir_global(0) * growth_increment,
                             dir_global(1) * growth_increment,
                             dir_global(2) * growth_increment);
    _active_front_node_growth_vectors.push_back(
        std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
  }
}

Real
MeshCut2DCCGUserObject::getTotalCrackLength() const
{
  return _total_crack_length;
}
