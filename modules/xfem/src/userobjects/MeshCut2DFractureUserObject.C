//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshCut2DFractureUserObject.h"

#include "XFEMFuncs.h"
#include "MooseError.h"
#include "MooseMesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_tools.h"
#include "algorithm"

#include "CrackFrontDefinition.h"

registerMooseObject("XFEMApp", MeshCut2DFractureUserObject);

InputParameters
MeshCut2DFractureUserObject::validParams()
{
  InputParameters params = MeshCut2DUserObjectBase::validParams();
  params.addClassDescription("XFEM mesh cutter for 2D models that defines cuts with a"
                             "mesh and uses fracture integrals to determine growth");
  params.addRequiredParam<Real>("growth_increment",
                                "Length to grow crack if k>k_critical or stress>stress_threshold");
  params.addParam<Real>("k_critical", "Critical fracture toughness.");
  params.addParam<Real>("stress_threshold", "Stress threshold for growing crack");
  params.addParam<VectorPostprocessorName>(
      "ki_vectorpostprocessor", "II_KI_1", "The name of the vectorpostprocessor that contains KI");
  params.addParam<VectorPostprocessorName>("kii_vectorpostprocessor",
                                           "II_KII_1",
                                           "The name of the vectorpostprocessor that contains KII");
  params.addParam<VectorPostprocessorName>(
      "stress_vectorpostprocessor",
      "The name of the vectorpostprocessor that contains crack front stress");
  params.addParam<std::string>("stress_vector_name",
                               "crack_tip_stress",
                               "The name of the stress vector in the stress_vectorpostprocessor");
  params.addParam<VectorPostprocessorName>(
      "k_critical_vectorpostprocessor",
      "The name of the vectorpostprocessor that contains critical fracture toughness at crack tip");
  params.addParam<std::string>(
      "k_critical_vector_name",
      "The name of the k_critical vector in the k_critical_vectorpostprocessor");
  return params;
}

MeshCut2DFractureUserObject::MeshCut2DFractureUserObject(const InputParameters & parameters)
  : MeshCut2DUserObjectBase(parameters),
    _growth_increment(getParam<Real>("growth_increment")),
    _use_k(isParamValid("k_critical") || isParamValid("k_critical_vectorpostprocessor")),
    _use_stress(isParamValid("stress_threshold")),
    _k_critical(isParamValid("k_critical") ? getParam<Real>("k_critical")
                                           : std::numeric_limits<Real>::max()),
    _stress_threshold(_use_stress ? getParam<Real>("stress_threshold")
                                  : std::numeric_limits<Real>::max()),
    _ki_vpp(_use_k ? &getVectorPostprocessorValue(
                         "ki_vectorpostprocessor",
                         getParam<VectorPostprocessorName>("ki_vectorpostprocessor"))
                   : nullptr),
    _kii_vpp(_use_k ? &getVectorPostprocessorValue(
                          "kii_vectorpostprocessor",
                          getParam<VectorPostprocessorName>("kii_vectorpostprocessor"))
                    : nullptr),
    _stress_vpp(_use_stress
                    ? &getVectorPostprocessorValue("stress_vectorpostprocessor",
                                                   getParam<std::string>("stress_vector_name"))
                    : nullptr),
    _k_critical_vpp(
        isParamValid("k_critical_vectorpostprocessor")
            ? &getVectorPostprocessorValue("k_critical_vectorpostprocessor",
                                           getParam<std::string>("k_critical_vector_name"))
            : nullptr)
{
  if (!_use_k && !_use_stress)
    paramError("k_critical",
               "Must set crack extension criterion with k_critical, k_critical_vectorpostprocessor "
               "or stress_threshold.");

  if (isParamValid("k_critical") && isParamValid("k_critical_vectorpostprocessor"))
    paramError("k_critical",
               "Fracture toughness cannot be specified by both k_critical and "
               "k_critical_vectorpostprocessor.");
}

void
MeshCut2DFractureUserObject::initialize()
{
  std::cout << "********MeshCut2DFractureUserObject::initialize()\n";
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
MeshCut2DFractureUserObject::isCutterMeshChanged() const
{
  return _is_mesh_modified;
}

void
MeshCut2DFractureUserObject::findActiveBoundaryGrowth()
{
  if ((!_ki_vpp || _ki_vpp->size() == 0) && (!_stress_vpp || _stress_vpp->size() == 0))
    return;

  if (_use_k && ((_ki_vpp->size() != _kii_vpp->size()) ||
                 (_ki_vpp->size() != _original_and_current_front_node_ids.size())))
    mooseError("ki_vectorpostprocessor and kii_vectorpostprocessor should have the same number of "
               "crack tips as CrackFrontDefinition.",
               "\n  ki size = ",
               _ki_vpp->size(),
               "\n  kii size = ",
               _kii_vpp->size(),
               "\n  cracktips in MeshCut2DFractureUserObject = ",
               _original_and_current_front_node_ids.size());

  if (_use_stress && ((_stress_vpp->size() != _original_and_current_front_node_ids.size())))
    mooseError("stress_vectorpostprocessor should have the same number of crack front points as "
               "CrackFrontDefinition.",
               "\n  stress_vectorpostprocessor size = ",
               _stress_vpp->size(),
               "\n  cracktips in MeshCut2DFractureUserObject = ",
               _original_and_current_front_node_ids.size());

  if (_k_critical_vpp && (_k_critical_vpp->size() != _original_and_current_front_node_ids.size()))
    mooseError("k_critical_vectorpostprocessor must have the same number of crack front points as "
               "CrackFrontDefinition.",
               "\n  k_critical_vectorpostprocessor size = ",
               _k_critical_vpp->size(),
               "\n  cracktips in MeshCut2DFractureUserObject = ",
               _original_and_current_front_node_ids.size());

  _active_front_node_growth_vectors.clear();

  for (unsigned int i = 0; i < _original_and_current_front_node_ids.size(); ++i)
  {
    bool was_crack_extended_kcrit = false;

    if (_use_k)
    {
      Real k_crit = _k_critical;
      if (_k_critical_vpp)
        k_crit = std::min(_k_critical_vpp->at(i), _k_critical);

      const Real ki = _ki_vpp->at(i);
      const Real kii = _kii_vpp->at(i);
      const Real k_sq = ki * ki + kii * kii;

      if (k_sq > k_crit * k_crit)
      {
        was_crack_extended_kcrit = true;

        const Real k_norm = std::sqrt(k_sq);
        const Real eps = std::max(1e-14, 1e-12 * k_norm);

        Real theta = 0.0;

        if (k_norm < eps)
        {
          theta = 0.0;
        }
        else if (std::abs(kii) < eps)
        {
          // Near pure Mode I:
          // KI > 0 -> forward growth, KI < 0 -> suppress or set theta = pi depending on model
          if (ki >= 0.0)
            theta = 0.0;
          else
            was_crack_extended_kcrit = false;
        }
        else
        {
          const Real root = std::sqrt(ki * ki + 8.0 * kii * kii);

          const Real theta_m = 2.0 * std::atan((ki - root) / (4.0 * kii));
          const Real theta_p = 2.0 * std::atan((ki + root) / (4.0 * kii));

          const Real sigma_tt_m =
              ki * (3.0 * std::cos(theta_m / 2.0) + std::cos(3.0 * theta_m / 2.0)) +
              kii * (-3.0 * std::sin(theta_m / 2.0) - 3.0 * std::sin(3.0 * theta_m / 2.0));

          const Real sigma_tt_p =
              ki * (3.0 * std::cos(theta_p / 2.0) + std::cos(3.0 * theta_p / 2.0)) +
              kii * (-3.0 * std::sin(theta_p / 2.0) - 3.0 * std::sin(3.0 * theta_p / 2.0));

          theta = (sigma_tt_m > sigma_tt_p) ? theta_m : theta_p;
        }

        if (was_crack_extended_kcrit)
        {
          RealVectorValue dir_cfc;
          dir_cfc(0) = std::cos(theta);
          dir_cfc(1) = std::sin(theta);
          dir_cfc(2) = 0.0;

          RealVectorValue dir_global =
              _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);

          const Real norm_dir = dir_global.norm();
          if (norm_dir > libMesh::TOLERANCE)
            dir_global /= norm_dir;
          else
            mooseError("Computed crack growth direction has near-zero norm.");

          Point nodal_offset(dir_global(0), dir_global(1), dir_global(2));
          nodal_offset *= _growth_increment;

          _active_front_node_growth_vectors.push_back(
              std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
        }
      }
    }

    if (_use_stress && !was_crack_extended_kcrit && _stress_vpp->at(i) > _stress_threshold)
    {
      RealVectorValue dir_cfc(1.0, 0.0, 0.0);
      RealVectorValue dir_global =
          _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);

      const Real norm_dir = dir_global.norm();
      if (norm_dir > libMesh::TOLERANCE)
        dir_global /= norm_dir;
      else
        mooseError("Computed crack growth direction has near-zero norm.");

      Point nodal_offset(dir_global(0), dir_global(1), dir_global(2));
      nodal_offset *= _growth_increment;

      _active_front_node_growth_vectors.push_back(
          std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
    }
  }
}

// void
// MeshCut2DFractureUserObject::findActiveBoundaryGrowth()
// {
//   // The k*_vpp & stress_vpp are empty (but not a nullptr) on the very first time step because this
//   // UO is called before the InteractionIntegral or crackFrontStress vpp
//   if ((!_ki_vpp || _ki_vpp->size() == 0) && (!_stress_vpp || _stress_vpp->size() == 0))
//     return;

//   if (_use_k && ((_ki_vpp->size() != _kii_vpp->size()) ||
//                  (_ki_vpp->size() != _original_and_current_front_node_ids.size())))
//     mooseError("ki_vectorpostprocessor and kii_vectorpostprocessor should have the same number of "
//                "crack tips as CrackFrontDefinition.",
//                "\n  ki size = ",
//                _ki_vpp->size(),
//                "\n  kii size = ",
//                _kii_vpp->size(),
//                "\n  cracktips in MeshCut2DFractureUserObject = ",
//                _original_and_current_front_node_ids.size());

//   if (_use_stress && ((_stress_vpp->size() != _original_and_current_front_node_ids.size())))
//     mooseError("stress_vectorpostprocessor should have the same number of crack front points as "
//                "CrackFrontDefinition.",
//                "\n  stress_vectorpostprocessor size = ",
//                _stress_vpp->size(),
//                "\n  cracktips in MeshCut2DFractureUserObject = ",
//                _original_and_current_front_node_ids.size());

//   if (_k_critical_vpp && ((_k_critical_vpp->size() != _original_and_current_front_node_ids.size())))
//     mooseError("k_critical_vectorpostprocessor must have the same number of crack front points as "
//                "CrackFrontDefinition.",
//                "\n  k_critical_vectorpostprocessor size = ",
//                _k_critical_vpp->size(),
//                "\n  cracktips in MeshCut2DFractureUserObject = ",
//                _original_and_current_front_node_ids.size());

//   _active_front_node_growth_vectors.clear();
//   for (unsigned int i = 0; i < _original_and_current_front_node_ids.size(); ++i)
//   {
//     // only extend crack with kcrit or nonlocal stress, never both.
//     bool was_crack_extended_kcrit = false;
//     if (_use_k)
//     {
//       Real k_crit = _k_critical;
//       if (_k_critical_vpp)
//         k_crit = std::min(_k_critical_vpp->at(i), _k_critical);

//       Real k_squared = _ki_vpp->at(i) * _ki_vpp->at(i) + _kii_vpp->at(i) * _kii_vpp->at(i);
//         std::cout<< "*** findActiveBoundaryGrowth(), debug 1, " << k_squared << ", " <<
//         _ki_vpp->at(i) << ", "<<  _kii_vpp->at(i)<< "," <<k_crit << "\n";
//       if (k_squared > (k_crit * k_crit) && _ki_vpp->at(i) > 0)
//       {
//         //     std::cout<< "*** findActiveBoundaryGrowth(), debug 2\n";
//         was_crack_extended_kcrit = true;
//         // growth direction in crack front coord (cfc) system based on the  max hoop stress
//         // criterion// criterion (maximum hoop stress)
//         // robust against KII -> 0 and avoids singular division in your current formula
//         Real ki = _ki_vpp->at(i);
//         Real kii = _kii_vpp->at(i);

//         Real theta = 0.0;

//         // scale-aware tolerance (tune if you like)
//         const Real k_norm = std::sqrt(ki * ki + kii * kii);
//         const Real eps = std::max(1e-14, 1e-12 * k_norm);

//         if (std::abs(kii) < eps || k_norm < eps)
//         {
//           // Pure Mode I (or essentially no driving): go straight
//           theta = 0.0;
//         }
//         else
//         {
//           // Use stable branch selection via sign(KII)
//           const Real sgn = (kii >= 0.0) ? 1 : -1;
//           const Real r = ki / kii; // safe because |kii| >= eps
//           Real sqrt_k = std::sqrt(ki * ki + 8*kii * kii);


//           //const Real a = 0.25 * (r - root); // stable expression
//           theta =  2 * std::atan((ki + sgn*sqrt_k) / (4 * kii));
//         }

//         // direction from angle
//         RealVectorValue dir_cfc;
//         dir_cfc(0) = std::cos(theta);
//         dir_cfc(1) = std::sin(theta);
//         dir_cfc(2) = 0.0;

//         // growth direction in global coord system based on the max hoop stress criterion
//         RealVectorValue dir_global =
//             _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);
//         Point dir_global_pt(dir_global(0), dir_global(1), dir_global(2));
//         Point nodal_offset = dir_global_pt * _growth_increment;
//         _active_front_node_growth_vectors.push_back(
//             std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
//       }
//     }
//     if (_use_stress && !was_crack_extended_kcrit && _stress_vpp->at(i) > _stress_threshold)
//     {
//       // crack will only be extended if it was not already extended by kcrit
//       // just extending the crack in the same direction it was going
//       RealVectorValue dir_cfc(1.0, 0.0, 0.0);
//       RealVectorValue dir_global =
//           _crack_front_definition->rotateFromCrackFrontCoordsToGlobal(dir_cfc, i);
//       Point dir_global_pt(dir_global(0), dir_global(1), dir_global(2));
//       Point nodal_offset = dir_global_pt * _growth_increment;
//       _active_front_node_growth_vectors.push_back(
//           std::make_pair(_original_and_current_front_node_ids[i].second, nodal_offset));
//     }
//   }
// }
