//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshCut2DCCGUserObject.h"

#include "XFEMFuncs.h"
#include "MooseError.h"
#include "MooseMesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_tools.h"

#include "CrackFrontDefinition.h"

registerMooseObject("XFEMApp", MeshCut2DCCGUserObject);

InputParameters
MeshCut2DCCGUserObject::validParams()
{
  InputParameters params = MeshCut2DUserObjectBase::validParams();
  params.addClassDescription("XFEM mesh cutter for 2D models that defines cuts with a"
                             "mesh and uses fracture integrals to determine growth");
  params.addRequiredParam<Real>("growth_increment",
                                "Length to grow crack if k>k_critical or stress>stress_threshold");
  params.addParam<Real>("paris_coeff", "Heading coefficient in Paris's law.");
  params.addParam<Real>("paris_exponent", "Exponential coefficient in Paris's law.");
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
    _growth_increment(getParam<Real>("growth_increment")),
    _ki_vpp(&getVectorPostprocessorValue(
        "ki_vectorpostprocessor", getParam<VectorPostprocessorName>("ki_vectorpostprocessor"))),
    _kii_vpp(&getVectorPostprocessorValue(
        "kii_vectorpostprocessor", getParam<VectorPostprocessorName>("kii_vectorpostprocessor"))),
    _c_vpp(getVectorPostprocessorValue("c_vectorpostprocessor",
                                       getParam<VectorPostprocessorName>("c_vectorpostprocessor"))),
    _total_crack_length(0),
    _paris_coeff(getParam<Real>("paris_coeff")),
    _paris_exponent(getParam<Real>("paris_exponent"))
{
  ci_old = 0.0;
  std::cout<<"_c_vpp size: " << _c_vpp.size() <<"\n";
}

void
MeshCut2DCCGUserObject::initialize()
{
  _is_mesh_modified = false;
  findActiveBoundaryGrowth();
  growFront();
  addNucleatedCracksToMesh();
  // update _crack_front_definition with nucleated nodes
  _crack_front_definition->updateNumberOfCrackFrontPoints(
      _original_and_current_front_node_ids.size());
  _crack_front_definition->isCutterModified(_is_mesh_modified);
  if (_is_mesh_modified)
    _crack_front_definition->updateCrackFrontPoints();

  std::vector<Point> p = getCrackFrontPoints(1);
  std::cout << "********MeshCut2DFractureUserObject::initialize(), crack tip location: " << p[0]
            << "\n";
}
bool
MeshCut2DCCGUserObject::isCutterMeshChanged() const
{
  return _is_mesh_modified;
}

void
MeshCut2DCCGUserObject::findActiveBoundaryGrowth()
{
  // k1 is empty (but not a nullptr) on the very first time step because this UO is called before
  // the InteractionIntegral or crackFrontStress vpp
  if ((_c_vpp.size() == 0))
    return;

  if (((_ki_vpp->size() != _kii_vpp->size()) ||
       (_ki_vpp->size() != _original_and_current_front_node_ids.size())))
    mooseError("ki_vectorpostprocessor and kii_vectorpostprocessor should have the same number of "
               "crack tips as CrackFrontDefinition.",
               "\n  ki size = ",
               _ki_vpp->size(),
               "\n  kii size = ",
               _kii_vpp->size(),
               "\n  cracktips in MeshCut2DFractureUserObject = ",
               _original_and_current_front_node_ids.size());

  _active_front_node_growth_vectors.clear();
  // std::cout<< "MeshCut2D *******************\n";
  // std::cout<< "MeshCut2DCCGUserOBject **********************, dt: " << _dt << "\n";
  double factor = _paris_exponent / (_paris_exponent + 1.0);
  for (unsigned int i = 0; i < _original_and_current_front_node_ids.size(); ++i)
  {

    if ( _c_vpp.at(i) > 0 )//&&
    //    _c_vpp.size() > 0) //(_c_vpp.size() > 0 && _c_vpp.at(i) > 0) //
    {
      Real ci = _c_vpp.at(i);
      Real growth_increment = _growth_increment;
      std::cout << "----update ci: " << ci << "\n";
      growth_increment =_paris_coeff * std::pow(ci, _paris_exponent) * _dt;

      _total_crack_length += growth_increment;
      std::cout << "growth_increment: " << growth_increment << "\n";

      // else
      //     growth_increment = std::pow(ci / 1e6 / 3600, 0.7) * 20 * 3600 * _dt;
      // if (_t > 2.0)
      //   growth_increment = std::pow(ci / 1e6 / 3600, 0.7) * 20 * 3600 * _dt;
      // growth direction in crack front coord (cfc) system based on the  max hoop stress
      // criterion
      Real ki = _ki_vpp->at(i);
      Real kii = _kii_vpp->at(i);
      Real sqrt_k = std::sqrt(ki * ki + 8.0 * kii * kii);
      Real theta = 2 * std::atan((ki - sqrt_k) / (4 * kii));
      RealVectorValue dir_cfc;
      dir_cfc(0) = std::cos(theta);
      dir_cfc(1) = std::sin(theta);
      dir_cfc(2) = 0;

      // growth direction in global coord system based on the max hoop stress criterion
      RealVectorValue dir_global =
          _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);
      Point dir_global_pt(dir_global(0), dir_global(1), dir_global(2));
      Point nodal_offset = dir_global_pt * growth_increment;
      _active_front_node_growth_vectors.push_back(
          std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
    }
  }
}

Real
MeshCut2DCCGUserObject::get_total_crack_length() const
{
  return _total_crack_length;
}
