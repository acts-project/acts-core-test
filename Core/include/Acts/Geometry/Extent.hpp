// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>
#include <vector>
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

using range_type = std::pair<double, double>;

// @brief Extent in space
///
/// This is a nested struct to the GeometryObject representation
/// which can be retrieved and used for surface parsing and will
/// give you the maximal extent in 3D space/
struct Extent {
  /// Possible maximal value
  static constexpr double maxval = std::numeric_limits<double>::max();

  /// Start value
  static constexpr range_type maxrange = {maxval, -maxval};

  // The different ranges
  std::vector<range_type> ranges = std::vector<range_type>(9, maxrange);

  /// Check the vertex
  /// @param vtx the Vertex to be checked
  void check(const Vector3D& vtx) {
    // min/max value check
    auto minMax = [&](BinningValue bval, double value) -> void {
      ranges[bval].first = std::min(value, ranges[bval].first);
      ranges[bval].second = std::max(value, ranges[bval].second);
    };
    // Walk through the binnin parameters
    for (int bval = 0; bval < binValues; ++bval) {
      BinningValue bValue = (BinningValue)bval;
      minMax(bValue, VectorHelpers::cast(vtx, bValue));
    }
  }
};
}  // namespace Acts