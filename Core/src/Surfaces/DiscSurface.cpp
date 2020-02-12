// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscSurface.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::DiscSurface::DiscSurface(const DiscSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(const GeometryContext& gctx,
                               const DiscSurface& other,
                               const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double rmin, double rmax, double hphisec)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const RadialBounds>(rmin, rmax, hphisec)) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double minhalfx, double maxhalfx, double maxR,
                               double minR, double avephi, double stereo)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const DiscTrapezoidBounds>(
          minhalfx, maxhalfx, maxR, minR, avephi, stereo)) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const DiscBounds> dbounds)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::move(dbounds)) {}

Acts::DiscSurface::DiscSurface(const std::shared_ptr<const DiscBounds>& dbounds,
                               const DetectorElementBase& detelement)
    : GeometryObject(), Surface(detelement), m_bounds(dbounds) {
  throw_assert(dbounds, "nullptr as DiscBounds");
}

Acts::DiscSurface& Acts::DiscSurface::operator=(const DiscSurface& other) {
  if (this != &other) {
    Acts::Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::Surface::SurfaceType Acts::DiscSurface::type() const {
  return Surface::Disc;
}

void Acts::DiscSurface::localToGlobal(const GeometryContext& gctx,
                                      const Vector2D& lposition,
                                      const Vector3D& /*gmom*/,
                                      Vector3D& position) const {
  // create the position in the local 3d frame
  Vector3D loc3Dframe(lposition[Acts::eLOC_R] * cos(lposition[Acts::eLOC_PHI]),
                      lposition[Acts::eLOC_R] * sin(lposition[Acts::eLOC_PHI]),
                      0.);
  // transport it to the globalframe (very unlikely that this is not needed)
  position = transform(gctx) * loc3Dframe;
}

bool Acts::DiscSurface::globalToLocal(const GeometryContext& gctx,
                                      const Vector3D& position,
                                      const Vector3D& /*gmom*/,
                                      Vector2D& lposition) const {
  // transport it to the globalframe (very unlikely that this is not needed)
  Vector3D loc3Dframe = (transform(gctx).inverse()) * position;
  lposition = Acts::Vector2D(perp(loc3Dframe), phi(loc3Dframe));
  return ((std::abs(loc3Dframe.z()) > s_onSurfaceTolerance) ? false : true);
}

const Acts::Vector2D Acts::DiscSurface::localPolarToLocalCartesian(
    const Vector2D& locpol) const {
  const DiscTrapezoidBounds* dtbo =
      dynamic_cast<const Acts::DiscTrapezoidBounds*>(&(bounds()));
  if (dtbo != nullptr) {
    double rMedium = dtbo->rCenter();
    double phi = dtbo->averagePhi();

    Vector2D polarCenter(rMedium, phi);
    Vector2D cartCenter = localPolarToCartesian(polarCenter);
    Vector2D cartPos = localPolarToCartesian(locpol);
    Vector2D Pos = cartPos - cartCenter;

    Acts::Vector2D locPos(
        Pos[Acts::eLOC_X] * sin(phi) - Pos[Acts::eLOC_Y] * cos(phi),
        Pos[Acts::eLOC_Y] * sin(phi) + Pos[Acts::eLOC_X] * cos(phi));
    return Vector2D(locPos[Acts::eLOC_X], locPos[Acts::eLOC_Y]);
  }
  return Vector2D(locpol[Acts::eLOC_R] * cos(locpol[Acts::eLOC_PHI]),
                  locpol[Acts::eLOC_R] * sin(locpol[Acts::eLOC_PHI]));
}

const Acts::Vector3D Acts::DiscSurface::localCartesianToGlobal(
    const GeometryContext& gctx, const Vector2D& lposition) const {
  Vector3D loc3Dframe(lposition[Acts::eLOC_X], lposition[Acts::eLOC_Y], 0.);
  return Vector3D(transform(gctx) * loc3Dframe);
}

const Acts::Vector2D Acts::DiscSurface::globalToLocalCartesian(
    const GeometryContext& gctx, const Vector3D& position,
    double /*unused*/) const {
  Vector3D loc3Dframe = (transform(gctx).inverse()) * position;
  return Vector2D(loc3Dframe.x(), loc3Dframe.y());
}

std::string Acts::DiscSurface::name() const {
  return "Acts::DiscSurface";
}

std::shared_ptr<Acts::DiscSurface> Acts::DiscSurface::clone(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return std::shared_ptr<DiscSurface>(this->clone_impl(gctx, shift));
}

Acts::DiscSurface* Acts::DiscSurface::clone_impl(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return new DiscSurface(gctx, *this, shift);
}

const Acts::SurfaceBounds& Acts::DiscSurface::bounds() const {
  if (m_bounds) {
    return (*(m_bounds.get()));
  }
  return s_noBounds;
}

Acts::PolyhedronRepresentation Acts::DiscSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg, bool triangulate) const {
  // Prepare vertices and faces
  std::vector<Vector3D> vertices;
  std::vector<std::vector<size_t>> faces;
  // Understand the disc
  bool fullDisc = m_bounds->coversFullAzimuth();
  bool toCenter = m_bounds->rMin() < s_onSurfaceTolerance;
  // If you have bounds you can create a polyhedron representation
  if (m_bounds) {
    // Helper function to create a face
    auto createFace = [&](std::vector<Vector2D>& vertices2D) -> void {
      // Split the vertices if you have
      std::vector<Vector3D> vertices3D;
      vertices3D.reserve(vertices2D.size());
      for (const auto& v : vertices2D) {
        vertices3D.push_back(transform(gctx) * Vector3D(v.x(), v.y(), 0.));
      }
      std::vector<size_t> face(vertices3D.size());
      std::iota(face.begin(), face.end(), vertices.size());
      vertices.insert(vertices.end(), vertices3D.begin(), vertices3D.end());
      faces.push_back(face);
    };
    auto v2D = m_bounds->vertices(lseg);
    // If not triangulate
    if (not triangulate) {
      // - just fill with 0 to N, split into two faces for rings
      if (fullDisc and not toCenter) {
        auto vsize = v2D.size();
        std::vector<Vector2D> half1(v2D.begin(), v2D.begin() + vsize / 4 + 1);
        half1.insert(half1.end(), v2D.begin() + 0.75 * vsize, v2D.end());
        half1.push_back(v2D[0.5 * vsize]);
        std::vector<Vector2D> half2(v2D.begin() + vsize / 4,
                                    v2D.begin() + vsize / 2);
        half2.push_back(v2D[0]);
        half2.insert(half2.end(), v2D.begin() + 0.5 * vsize,
                     v2D.begin() + 0.75 * vsize + 1);
        createFace(half1);
        createFace(half2);
      } else {
        createFace(v2D);
      }
    } else {
      /// Triangulation
      // Disc to the center are handled with a central anker point
      if (toCenter) {
        vertices.reserve(v2D.size() + 1);
        if (fullDisc) {
          vertices.push_back(Vector3D(0., 0., 0.));
        }
        for (const auto& v : v2D) {
          vertices.push_back(transform(gctx) * Vector3D(v.x(), v.y(), 0.));
          if (vertices.size() > 2) {
            std::vector<size_t> face = {0, vertices.size() - 2,
                                        vertices.size() - 1};
            faces.push_back(face);
          }
        }
        if (fullDisc) {
          faces.push_back({0, vertices.size() - 1, 1});
        }
      } else {
        // Reorder the vertices for triangulation
        vertices.reserve(v2D.size());
        for (size_t iv = 0; iv < 0.5 * v2D.size(); ++iv) {
          auto vf = v2D[iv];
          auto vs = v2D[v2D.size() - 1 - iv];
          vertices.push_back(transform(gctx) * Vector3D(vf.x(), vf.y(), 0.));
          vertices.push_back(transform(gctx) * Vector3D(vs.x(), vs.y(), 0.));
        }
        for (size_t iv = 0; iv < vertices.size() + 1; ++iv) {
          if (iv > 2) {
            std::vector<size_t> face = {iv - 3, iv - 2, iv - 1};
            if (iv % 2) {
              std::reverse(face.begin(), face.end());
            }
            faces.push_back(face);
          }
        }
        // Close if it's a ring
        if (fullDisc) {
          faces.push_back({0, 1, vertices.size() - 2});
          faces.push_back({vertices.size() - 2, 1, vertices.size() - 1});
        }
      }
    }
  } else {
    throw std::domain_error(
        "Polyhedron repr of boundless surface not possible.");
  }
  return PolyhedronRepresentation(vertices, faces);
}
