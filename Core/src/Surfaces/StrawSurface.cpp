// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Surfaces/detail/VertexHelper.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Surfaces/InfiniteBounds.hpp"

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 double radius, double halez)
    : GeometryObject(), LineSurface(std::move(htrans), radius, halez) {}

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 std::shared_ptr<const LineBounds> lbounds)
    : GeometryObject(), LineSurface(std::move(htrans), std::move(lbounds)) {}

Acts::StrawSurface::StrawSurface(
    const std::shared_ptr<const LineBounds>& lbounds,
    const DetectorElementBase& detelement)
    : GeometryObject(), LineSurface(lbounds, detelement) {}

Acts::StrawSurface::StrawSurface(const Acts::StrawSurface& other)
    : GeometryObject(), LineSurface(other) {}

Acts::StrawSurface::StrawSurface(const GeometryContext& gctx,
                                 const StrawSurface& other,
                                 const Transform3D& transf)
    : GeometryObject(), LineSurface(gctx, other, transf) {}

Acts::StrawSurface& Acts::StrawSurface::operator=(const StrawSurface& other) {
  if (this != &other) {
    LineSurface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

std::shared_ptr<Acts::StrawSurface> Acts::StrawSurface::clone(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return std::shared_ptr<StrawSurface>(this->clone_impl(gctx, shift));
}

Acts::StrawSurface* Acts::StrawSurface::clone_impl(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return new StrawSurface(gctx, *this, shift);
}

Acts::PolyhedronRepresentation Acts::StrawSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg, bool triangulate) const {
  std::vector<Vector3D> vertices;
  std::vector<std::vector<size_t>> faces;

  const Transform3D& ctransform = transform(gctx);
  double hZ = m_bounds->halflengthZ();

  // Draw the bounds if more than one segment are chosen
  if (lseg > 1) {
    double r = m_bounds->r();
    auto phiSegs = detail::VertexHelper::phiSegments();
    // Write the two bows/circles on either side
    std::vector<int> sides = {-1, 1};
    for (auto& side : sides) {
      for (size_t iseg = 0; iseg < phiSegs.size() - 1; ++iseg) {
        int addon = (iseg == phiSegs.size() - 2) ? 1 : 0;
        /// Helper method to create the segment
        detail::VertexHelper::createSegment(
            vertices, r, phiSegs[iseg], phiSegs[iseg + 1], lseg, addon,
            Vector3D(0., 0., side * hZ), ctransform);
      }
    }
    // Write the faces from the built vertices
    size_t nqfaces = 0.5 * (vertices.size() - 2);
    for (size_t iface = 0; iface < nqfaces; ++iface) {
      size_t p2 = (iface + 1 == nqfaces) ? 0 : iface + 1;
      if (not triangulate) {
        std::vector<size_t> face = {iface, p2, p2 + nqfaces, nqfaces + iface};
        faces.push_back(face);
      } else {
        // Create two if configured to triangulate
        std::vector<size_t> triA = {iface, p2, p2 + nqfaces};
        faces.push_back(triA);
        std::vector<size_t> triB = {p2 + nqfaces, nqfaces + iface, iface};
        faces.push_back(triB);
      }
    }
  }

  size_t bvertices = vertices.size();
  Vector3D left(0, 0, -hZ);
  Vector3D right(0, 0, hZ);
  // The central wire/straw
  vertices.push_back(ctransform * left);
  vertices.push_back(ctransform * right);
  if (not triangulate) {
    faces.push_back({bvertices, bvertices + 1});
  } else {
    vertices.push_back(ctransform * Vector3D(0., 0., 0.));
    faces.push_back({bvertices, bvertices + 2, bvertices + 1});
  }

  return PolyhedronRepresentation(vertices, faces);
}
