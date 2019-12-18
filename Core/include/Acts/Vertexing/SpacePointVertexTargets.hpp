// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <utility>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief a phi/z target and cast with projection along rho
class CylindricalDetectorTarget {
 public:
  enum TargetType : int { vertex = -1, cylinder = 0, disk = 1 };

  /// Default contstructor for vertex target
  CylindricalDetectorTarget() : m_targetType(vertex) {}

  /// Constructor of a phi-z target
  ///
  /// @param target is either the radius or the z value (boolean steered)
  /// @param targetRange is the range this target covers
  /// @param layerID is the layer identifier for bit setting
  ///  required to be > 0 and unique
  CylindricalDetectorTarget(const Eigen::Array<double, 2, 2>& targetMatrix,
                            unsigned int layerID)
      : m_targetID(layerID) {
    auto rangeLengths =
        (targetMatrix.block<2, 1>(1, 0) - targetMatrix.block<2, 1>(0, 0));
    m_targetType = (rangeLengths(0, 0) > rangeLengths(0, 1) ? cylinder : disk);
    m_projectionType = (m_targetType == cylinder) ? 1 : 0;
    m_targetRange = targetMatrix.block<1, 2>(m_projectionType, 0);
    m_target = 0.5 * (m_targetRange(0, 0) + m_targetRange(1, 0));
    m_projectionRange = targetMatrix.block<1, 2>(m_targetType, 0);
    m_projectionScale =
        1. / (m_projectionRange(1, 0) - m_projectionRange(0, 0));
  }

  /// This casts the projection parameters
  ///
  /// @tparam input_t the templated input point
  ///
  template <typename input_t>
  auto projectionParameters(const input_t& point) const {
    Eigen::Array<double, 1, 2> projp;
    projp << point.rawValues(0, 0), point.rawValues(0, m_targetType + 1);
    return projp;
  }

  /// Get and fill the raw value types
  static auto rawValues(const Vector3D& position) {
    Eigen::Array<double, 1, 3> rValues;
    rValues << VectorHelpers::phi(position), position.z(),
        VectorHelpers::perp(position);
    return rValues;
  }

  /// This casts the projection parameters
  template <typename input_t>
  auto k(const input_t& point) const {
    return point.rawValues(0, m_projectionType + 1);
  }

  /// Double cast operator
  operator double() const { return m_target; }

  /// Scale parameters if required
  ///
  /// @param projp are the projected parameters at target (unscaled)
  void scaleProjection(Eigen::Array<double, 1, 2>& projp) const {
    if (m_targetType == vertex)
      return;
    // Cylinder scales to (-1,1), disk scales from (-1,0) or (0,1)
    // depending on target side
    projp(0, 1) =
        (-1. + m_targetType) + (projp(0, 1) - m_projectionRange(0, 0)) *
                                   (2 - m_targetType) * m_projectionScale;
    // Mirror the disk if necessary
    if (m_targetType == disk) {
      projp(0, 1) = std::copysign(projp(0, 1), m_target);
    }
  }

  /// The layer
  unsigned int targetID() const { return m_targetID; }

  /// Expose the type
  TargetType targetType() const { return m_targetType; }

 public:
  double m_target = 0.;                      //< The target radius
  Eigen::Array<double, 2, 1> m_targetRange;  //< The target range covered
  double m_projectionScale = 0.;             // The spread of the target reange
  Eigen::Array<double, 2, 1> m_projectionRange;  //< The target range covered
  double m_targetID = 0;                         //< The layer identification
  TargetType m_targetType = cylinder;
  unsigned int m_projectionType = 0;
};

}  // namespace Acts