// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Utilities/BinningType.hpp"

Acts::Polyhedron& Acts::Polyhedron::operator+=(const Acts::Polyhedron& other) {
  size_t cvert = vertices.size();
  vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());
  /// Add the new faces with offsets
  auto join = [&](std::vector<face_type>& existing,
                  const std::vector<face_type>& additional) -> void {
    for (const auto& aface : additional) {
      face_type nface = aface;
      std::transform(nface.begin(), nface.end(), nface.begin(),
                     [&](size_t x) { return (x + cvert); });
      existing.push_back(nface);
    }
  };
  // For faces and triangular mesh
  join(faces, other.faces);
  join(triangularMesh, other.triangularMesh);

  return (*this);
}

Acts::Extent Acts::Polyhedron::extent(const Transform3D& transform) const {
  Extent extent;
  for (const auto& vtx : vertices) {
    extent.check(transform * vtx);
  }
  // For purely planar surfaces we need to check if the (0,0,0) is
  // inside any of the convex faces
  if (std::abs(extent.ranges[binZ].first) < s_onSurfaceTolerance and
      std::abs(extent.ranges[binZ].second) < s_onSurfaceTolerance) {
    Vector3D origin = transform * Vector3D(0., 0., 0.);
    for (const auto& face : faces) {
      std::vector<Vector3D> tface;
      tface.reserve(face.size());
      for (auto f : face) {
        tface.push_back(transform * vertices[f]);
      }
      if (detail::VertexHelper::isInsidePolygon(origin, tface)) {
        extent.ranges[binR].first = 0.;
        break;
      }
    }
  }
  return extent;
}
