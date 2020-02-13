// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

Acts::Polyhedron& Acts::Polyhedron::operator+=(const Polyhedron& other);
{
  size_t cvert = vertices.size();
  vertices.insert(vertices.begin(), other.vertices.begin(),
                  other.vertices.end());
  std::cout << "Adding " << other.vertices.size() << " to " << cvert
            << std::endl;
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
}
