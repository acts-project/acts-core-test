// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CompoundLayer.hpp"

Acts::CompoundLayer::CompoundLayer(
    std::shared_ptr<Surface> repSurface, double repThickness,
    std::unique_ptr<CompoundSurfaceArrays> surfaceArrays,
    std::unique_ptr<CompoundApproachDescriptors> apporachDescriptors,
    LayerType laytyp)
    : Layer(nullptr, repThickness, nullptr, laytyp),
      m_representingSurface(std::move(repSurface)),
      m_surfaceArrays(std::move(surfaceArrays)),
      m_apporachDescriptors(std::move(apporachDescriptors)) {}

Acts::Polyhedron Acts::CompoundLayer::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg) const {
  Polyhedron phed;
  if (m_apporachDescriptors != nullptr) {
    for (const auto& aDescriptor : m_apporachDescriptors->arrayObjects()) {
      for (const auto& aSurface : aDescriptor->containedSurfaces()) {
        phed += aSurface->polyhedronRepresentation(gctx, lseg);
      }
    }
    return phed;
  }
  return surfaceRepresentation().polyhedronRepresentation(gctx, lseg);
}