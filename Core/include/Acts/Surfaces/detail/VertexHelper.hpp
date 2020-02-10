// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

namespace detail {

/// @brief Helpers to draw segments between phi min/max
namespace VertexHelper {

/// A method that inserts the cartesian extrema points and segments
/// a curved segment into sub segments
///
/// @param phiMin the minimum Phi of the bounds object
/// @param phiMax the maximum Phi of the bounds object
/// @param phiRef is a vector of reference phi values to be included as well
/// @param phiTolerance is the tolerance for reference phi insertion
///
/// @return a vector
static std::vector<double> phiSegments(double phiMin = -M_PI,
                                       double phiMax = M_PI,
                                       std::vector<double> phiRefs = {},
                                       double phiTolerance = 0.) {
  // This is to ensure that the extrema are built regardless of number
  // of segments
  std::vector<double> phiSegments;
  std::vector<double> quarters = {-M_PI, -0.5 * M_PI, 0., 0.5 * M_PI, M_PI};
  // It does not cover the full azimuth
  if (phiMin != -M_PI or phiMax != M_PI) {
    phiSegments.push_back(phiMin);
    for (unsigned int iq = 1; iq < 4; ++iq) {
      if (phiMin < quarters[iq] and phiMax > quarters[iq]) {
        phiSegments.push_back(quarters[iq]);
      }
    }
    phiSegments.push_back(phiMax);
  } else {
    phiSegments = quarters;
  }
  // Insert the reference phis if
  if (not phiRefs.empty()) {
    for (const auto& phiRef : phiRefs) {
      // Trying to find the right patch
      auto match = std::find_if(
          phiSegments.begin(), phiSegments.end(), [&](double phiSeg) {
            return std::abs(phiSeg - phiRef) < phiTolerance;
          });
      if (match == phiSegments.end()) {
        phiSegments.push_back(phiRef);
      }
    }
    std::sort(phiSegments.begin(), phiSegments.end());
  }
  return phiSegments;
}

}  // namespace VertexHelper

}  // namespace detail

}  // namespace Acts