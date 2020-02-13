// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

// Helper
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

// The class to test
#include <fstream>
#include "Acts/Geometry/Polyhedron.hpp"

namespace Acts {

using namespace UnitLiterals;

using IdentifiedPolyderon = std::pair<std::string, Polyhedron>;

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(PolyhedronTest) {
  std::vector<Vector3D> tvertices = {Vector3D(-1, -1, 0.), Vector3D(1., -1, 0.),
                                     Vector3D(0., 1., 0.)};
  std::vector<std::vector<size_t>> tfaces = {{0, 1, 2}};

  Polyhedron triangle(tvertices, tfaces, tfaces);
  BOOST_CHECK(tvertices == triangle.vertices);
  BOOST_CHECK(tfaces == triangle.faces);
  BOOST_CHECK(tfaces == triangle.triangularMesh);

  std::ofstream tStream;
  tStream.open("PolyhedronTriangle.obj");
  ObjHelper objtH;
  triangle.draw(objtH);
  objtH.write(tStream);
  tStream.close();

  std::vector<Vector3D> rvertices = {Vector3D(-1, -2, 0.), Vector3D(1., -2, 0.),
                                     Vector3D(1., -1., 0.),
                                     Vector3D(-1., -1., 0.)};
  std::vector<std::vector<size_t>> rfaces = {{0, 1, 2, 3}};
  std::vector<std::vector<size_t>> rmesh = {{0, 1, 2}, {2, 3, 0}};
  Polyhedron rectangle(rvertices, rfaces, rmesh);
  BOOST_CHECK(rvertices == rectangle.vertices);
  BOOST_CHECK(rfaces == rectangle.faces);
  BOOST_CHECK(rmesh == rectangle.triangularMesh);

  std::ofstream rStream;
  rStream.open("PolyhedronRectangle.obj");
  ObjHelper objrH;
  rectangle.draw(objrH);
  objrH.write(rStream);
  rStream.close();

  // Now add them
  Polyhedron tr();
  // tr += triangle;
  // BOOST_CHECK(tr.vertices == rectangle.vertices);
  // BOOST_CHECK(tr.faces == rectangle.faces);
  // BOOST_CHECK(tr.triangularMesh == rectangle.triangularMesh);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts