// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Utilities/ThrowAssert.hpp"

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               double radius, double halez)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const LineBounds>(radius, halez)) {}

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const LineBounds> lbounds)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::move(lbounds)) {}

Acts::LineSurface::LineSurface(const std::shared_ptr<const LineBounds>& lbounds,
                               const DetectorElementBase& detelement)
    : GeometryObject(), Surface(detelement), m_bounds(lbounds) {
  throw_assert(lbounds, "LineBounds must not be nullptr");
}

Acts::LineSurface::LineSurface(const LineSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::LineSurface::LineSurface(const GeometryContext& gctx,
                               const LineSurface& other,
                               const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::LineSurface& Acts::LineSurface::operator=(const LineSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::PolyhedronRepresentation Acts::LineSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg, bool /*ignored*/) const {
  std::vector<Vector3D> vertices;
  std::vector<std::vector<size_t>> faces;

  if (lseg <= 1) {
    throw std::domain_error(
        "Polyhedron repr of cylinder with 1 div is undefined");
  }

  if (m_bounds != nullptr) {
    auto linetra = transform(gctx);
    auto linedir = linetra.matrix().block<3, 1>(0, 2);
    auto linectr = center(gctx);
    // line vertices and counding cylinder
    vertices.reserve(2 * lseg + 2);
    faces.reserve(lseg + 1);

    double hlZ = m_bounds->halflengthZ();
    double r = m_bounds->r();
    double phiStep = 2 * M_PI / lseg;

    // The vertices for the facettes
    for (unsigned int iseg = 0; iseg < lseg; ++iseg) {
      double phi = -M_PI + iseg * phiStep;
      double cphi = std::cos(phi);
      double sphi = std::sin(phi);
      vertices.push_back(linetra * Vector3D(r * cphi, r * sphi, -hlZ));
      vertices.push_back(linetra * Vector3D(r * cphi, r * sphi, hlZ));
    }
    for (size_t v = 0; v < vertices.size() - 2; v = v + 2) {
      faces.push_back({v, v + 1, v + 3, v + 2});
    }
    faces.push_back({vertices.size() - 2, vertices.size() - 1, 1, 0});

    // close the facettes
    faces.push_back({2 * (lseg - 1), 0, 1, 2 * lseg});

    // The two line vertices
    vertices.push_back(linectr - hlZ * linedir);
    vertices.push_back(linectr + hlZ * linedir);
    faces.push_back({vertices.size() - 2, vertices.size() - 1});
  }
  return PolyhedronRepresentation(vertices, faces);
}
