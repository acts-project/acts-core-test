// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iostream>
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @struct ProtoLayer
///
/// Encapsulates min/max boundaries that will be turned into a layer.
/// The struct allows this information to be obtained in a consistent
/// way, or be caller provided.

struct ProtoLayer {
 public:
  /// The extent of the ProtoLayer
  Extent environment;

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The vector of surfaces to consider
  ProtoLayer(const GeometryContext& gctx,
             const std::vector<const Surface*>& surfaces);

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The vector of surfaces to consider
  ProtoLayer(const GeometryContext& gctx,
             const std::vector<std::shared_ptr<const Surface>>& surfaces);

  // Defaulated empty constructor
  ProtoLayer() = default;

  std::ostream& toStream(std::ostream& sl) const;

 private:
  /// Helper method which performs the actual min/max calculation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The surfaces to build this protolayer out of
  void measure(const GeometryContext& gctx,
               const std::vector<const Surface*>& surfaces);
};
}  // namespace Acts