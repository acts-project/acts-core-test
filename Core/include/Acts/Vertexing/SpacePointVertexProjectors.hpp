// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <utility>
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @brief Linear projection of a seed to a target which will
/// return the projected values at the target
struct LinearProjector {
  /// @brief projection method
  ///
  /// @tparam scan_point_t The projection point
  /// @tparam target_t The projection target
  ///
  /// @todo change to parameter pack or range
  ///
  /// @return inwards and outwards projections
  template <typename scan_point_t, typename target_t>
  auto project(const scan_point_t& p0, const scan_point_t& p1,
               const target_t& target) const {
    auto p0array = target.projectionParameters(p0);
    auto p1array = target.projectionParameters(p1);
    double k = (target.k(p0) - double(target)) / (target.k(p1) - target.k(p0));
    Eigen::Array<double, 1, 2> karray;
    karray << k, k;
    Eigen::Array<double, 1, 2> proj = (p0array - karray * (p1array - p0array));
    target.scaleProjection(proj);
    return proj;
  }
};

}  // namespace Acts
