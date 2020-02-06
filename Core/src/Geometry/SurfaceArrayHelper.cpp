// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayHelper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Geometry/SurfaceArrayHelper.hpp"
#include <algorithm>
#include <cmath>
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::SortedSurfaceVector Acts::SurfaceArrayHelper::cylinders(
    const GeometryContext& gctx, const SurfaceVector& surfaces,
    const Options& options) const {
  // The split value
  BinningValue splitBin = binR;
  double splitTolerance = options.rSplitTolerance;
  // Return the sorted
  auto sSurfaces = sort(gctx, surfaces, splitBin, splitTolerance);
  // Estimate the binning

  return sSurfaces;
}

Acts::SortedSurfaceVector Acts::SurfaceArrayHelper::discs(
    const GeometryContext& gctx, const SurfaceVector& surfaces,
    const Options& options) const {
  // The split value
  BinningValue splitBin = binZ;
  double splitTolerance = options.zSplitTolerance;
  // Return the sorted
  auto sSurfaces = sort(gctx, surfaces, splitBin, splitTolerance);
  //

  // Estimate the binning
  return sSurfaces;
}

Acts::SortedSurfaceVector Acts::SurfaceArrayHelper::planes(
    const GeometryContext& gctx, const SurfaceVector& surfaces,
    const Options& options) const {
  // The split value for planes is Z (for the moment),
  // as local coordinate system is x/y
  BinningValue splitBin = binZ;
  double splitTolerance = options.zSplitTolerance;
  // Return the sorted
  auto sSurfaces = sort(gctx, surfaces, splitBin, splitTolerance);

  // Estimate the binning
  return sSurfaces;
}

Acts::SortedSurfaceVector Acts::SurfaceArrayHelper::sort(
    const GeometryContext& gctx, const SurfaceVector& surfaces,
    BinningValue splitBin, double splitTolerance) const {
  SortedSurfaceVector sSurfaces;
  for (auto& sf : surfaces) {
    double cvalue = sf->binningPositionValue(gctx, splitBin);
    // Trying to find the right patch
    auto patch = std::find_if(
        sSurfaces.begin(), sSurfaces.end(), [&](SurfacePatch& surfacePatch) {
          return std::abs(surfacePatch.splitValue - cvalue) < splitTolerance;
        });
    /// patch not found, assign
    if (patch == sSurfaces.end()) {
      // Create a new patch and fill it
      SurfacePatch newPatch;
      newPatch.splitValue = cvalue;
      newPatch.splitBin = splitBin;
      newPatch.splitTolerance = splitTolerance;
      // Take and assign
      sSurfaces.push_back(newPatch);
      patch = (sSurfaces.end() - 1);
    }
    patch->surfaces.push_back(sf);
  }
  // Sort the patches
  std::sort(sSurfaces.begin(), sSurfaces.end());

  return sSurfaces;
}
