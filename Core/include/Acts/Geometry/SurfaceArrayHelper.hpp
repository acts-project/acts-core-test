// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayHelper.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <climits>
#include <vector>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class Surface;

using SurfaceVector = std::vector<const Surface*>;
auto mval = std::numeric_limits<double>::max();

/// @brief helper return struct for dimensions, binning & surfaces
struct SurfacePatch {
  /// x min/max from parsing
  double xMin = mval, xMax = -mval;
  /// y min/max from parsing
  double yMin = mval, yMax = -mval;
  /// r min/max from parsing
  double rMin = mval, rMax = -mval;
  /// z min/max from parsing
  double zMin = mval, zMax = -mval;
  /// Reference binning for sorting
  BinningValue splitBin = binR;
  /// Reference binning value for sorting
  double splitValue = 0.;
  /// Split tolerance
  double splitTolerance = 0.;

  /// number of bins
  unsigned int bins0 = 1, bins1 = 1;

  /// Surfaces needed for the array
  SurfaceVector surfaces;

  /// <operator for sorting
  bool operator<(const SurfacePatch& other) const {
    return (splitValue < other.splitValue);
  }
};

using SortedSurfaceVector = std::vector<SurfacePatch>;

/// @brief Helper tool that sorts sesnsitive vectors into layer
/// sorted packages in order to build them into Layers
///
/// It also performs ring layout detection, i.e. detects if a
/// disk setup is better described by rings that are sorted into
/// individual volumes.
///
class SurfaceArrayHelper {
 public:
  /// Nested configuraiton struct
  struct Config {};

  /// Nested options struct
  struct Options {
    /// Tolerance to split in r
    double rSplitTolerance = 0.;
    /// Tolerance to split in z
    double zSplitTolerance = 0.;
  };

  /// Constructor with explicit config
  ///
  /// @param cfg Explicit config struct
  /// @param logger logging instance
  SurfaceArrayHelper(const Config& cfg,
                     std::unique_ptr<const Logger> logger =
                         getDefaultLogger("SurfaceArrayCreator", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  /// Destructor
  virtual ~SurfaceArrayHelper() = default;

  /// Sort surfaces into cylindrical patches
  /// @param gctx The geometry context
  /// @param surfaces The input surfaces (unordered)
  /// @param options The sorting options
  SortedSurfaceVector cylinders(const GeometryContext& gctx,
                                const SurfaceVector& surfaces,
                                const Options& options) const;

  /// Sort surfaces into cylindrical patches
  /// @param gctx The geometry context
  /// @param surfaces The input surfaces (unordered)
  /// @param options The sorting options
  SortedSurfaceVector discs(const GeometryContext& gctx,
                            const SurfaceVector& surfaces,
                            const Options& options) const;

  /// Sort surfaces into planar patches
  /// @param gctx The geometry context
  /// @param surfaces The input surfaces (unordered)
  /// @param options The sorting options
  SortedSurfaceVector planes(const GeometryContext& gctx,
                             const SurfaceVector& surfaces,
                             const Options& options) const;

  /// Set logging instance
  /// @param logger is the logging instance to be set
  void setLogger(std::unique_ptr<const Logger> logger) {
    m_logger = std::move(logger);
  }

 private:
  /// configuration object
  Config m_cfg;

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logger
  const Logger& logger() const { return *m_logger; }

  /// Parse and sort boundaries
  /// @param surfaces The input surfaces (unordered)
  /// @param splitBin The binning value for the split
  /// @param splitTolerance The double value tolerance to split
  SortedSurfaceVector sort(const GeometryContext& gctx,
                           const SurfaceVector& surfaces, BinningValue splitBin,
                           double splitTolerance) const;
};

}  // namespace Acts