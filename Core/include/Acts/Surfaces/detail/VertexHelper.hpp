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
                                       double phiTolerance = 1e-6) {
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

/// Helper method to create a regular 2 or 3 D segment
///  between two phi values
///
/// @tparam vertex_t Type of vertex to be applied
/// @tparam transform_t Optional transform
///
/// @param vertices [in,out] The 3D vertices to be filled
/// @param r The radius
/// @param phi1 The first phi value
/// @param phi2 The second phi value
/// @param lseg The number of segments for full 2*PI
/// @param addon The additional segments to be built
/// @param offset The out of plane offset position of the bow
/// @param transform The transform applied (optional)
template <typename vertex_t, typename transform_t>
void createSegment(std::vector<vertex_t>& vertices, double r, double phi1,
                   double phi2, unsigned int lseg, int addon = 0,
                   const vertex_t& offset = vertex_t::Zero(),
                   const transform_t& transform = transform_t::Identity()) {
  // Calculate the number of segments - 1 is the minimum
  unsigned int segs = std::abs(phi2 - phi1) / (2 * M_PI) * lseg;
  segs = segs > 0 ? segs : 1;
  double phistep = (phi2 - phi1) / segs;
  // Create the segments
  for (unsigned int iphi = 0; iphi < segs + addon; ++iphi) {
    double phi = phi1 + iphi * phistep;
    vertex_t vertex = vertex_t::Zero();
    vertex(0) = r * std::cos(phi);
    vertex(1) = r * std::sin(phi);
    vertex = vertex + offset;
    vertices.push_back(transform * vertex);
  }
};

}  // namespace VertexHelper

}  // namespace detail

}  // namespace Acts