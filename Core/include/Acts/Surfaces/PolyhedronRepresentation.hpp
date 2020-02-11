// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <vector>

#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @class PolyhedronRepresentation
///
/// Struct which contains a cartesian approximation for any surface type.
/// It contains a list of cartesian vertices in the global frame, and
/// additionally
/// a list of lists of indices which indicate which vertices form a face.
/// Each entry in @c faces is a face, which is in turn a list of vertices
/// that need to be connected to form a face.
/// This allows the @c objString method to produce a ready-to-go obj output.
struct PolyhedronRepresentation {
  /// Default constructor from a vector of vertices and a vector of faces
  /// @param verticesIn The 3D global vertices that make up the object
  /// @param facesIn List of lists of indices for faces.
  /// @note This creates copies of the input vectors
  PolyhedronRepresentation(const std::vector<Vector3D>& verticesIn,
                           const std::vector<std::vector<size_t>>& facesIn)
      : vertices(verticesIn), faces(facesIn) {}

  /// list of 3D vertices as vectors
  std::vector<Vector3D> vertices;

  /// list of faces connecting the vertices.
  /// each face is a list of vertices v
  /// corresponding to the vertex vector above
  std::vector<std::vector<size_t>> faces;

  /// Draw method for polyhedrons
  ///
  /// @tparam helper_t The draw helper
  ///
  /// @param helper The draw helper object (visitor pattern)
  /// @param decompose Boolean that forces a decomposition into
  /// individual vertices
  template <typename helper_t>
  void draw(helper_t& helper, bool decompose = false) const {
    // vertices and faces are
    if (not decompose) {
      helper.faces(vertices, faces);
    } else {
      for (const auto& face : faces) {
        std::vector<Vector3D> face_vtx;
        for (size_t i : face) {
          face_vtx.push_back(vertices[i]);
        }
        helper.face(face_vtx);
      }
    }
  }

  /// Maximum extent of the surface in space
  ///
  /// @param transform is the an (optional) transform
  /// to apply to the vertices for estimation the extent
  /// with respect to a coordinate frame
  ///
  /// @return ranges that describe the space taken by this surface
  GeometryObject::Extent surfaceExtent(
      const Transform3D& transform = Transform3D::Identity()) const {
    GeometryObject::Extent surfaceExtend;
    for (const auto& vtx : vertices) {
      surfaceExtend.check(transform * vtx);
    }
    return surfaceExtend;
  }
};
}  // namespace Acts
