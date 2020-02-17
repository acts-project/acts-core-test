// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <algorithm>
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace Acts {

ProtoLayer::ProtoLayer(const GeometryContext& gctx,
                       const std::vector<const Surface*>& surfaces) {
  measure(gctx, surfaces);
}

ProtoLayer::ProtoLayer(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<const Surface>>& surfaces) {
  measure(gctx, unpack_shared_vector(surfaces));
}

std::ostream& ProtoLayer::toStream(std::ostream& sl) const {
  sl << "ProtoLayer with dimensions (min/max)" << std::endl;
  sl << " - r : " << minR << " - " << envR.first << " / " << maxR << " + "
     << envR.second << std::endl;
  sl << " - z : " << minZ << " - " << envZ.first << " / " << maxZ << " + "
     << envZ.second << std::endl;
  sl << " - phi : " << minPhi << " - " << envPhi.first << " / " << maxPhi
     << " + " << envPhi.second << std::endl;

  return sl;
}

void ProtoLayer::measure(const GeometryContext& gctx,
                         const std::vector<const Surface*>& surfaces) {
  for (const auto& sf : surfaces) {
    // Take the thickness in account if necessary
    double thickness = 0;
    auto sfPolyhedron = sf->polyhedronRepresentation(gctx,1);
    const DetectorElementBase* element = sf->associatedDetectorElement();
    if (element != nullptr) {
      thickness = element->thickness();
      // We need a translation along and opposite half thickness
      Vector3D sfNormal = sf->normal(gctx, center(gctx));
      std::vector<Translation3D> deltaT = 
      { Translation3D(-0.5*thiskness*sfNormal),
        Translation3D(-0.5*thiskness*sfNormal)};
      for (const auto& dT : deltaT){
          environment += sfPolyhedron.extent(dT); 
      }  
      continue;
    } 
    environment += sfPolyhedron.extent(); 
  }                          
}

}  // namespace Acts
