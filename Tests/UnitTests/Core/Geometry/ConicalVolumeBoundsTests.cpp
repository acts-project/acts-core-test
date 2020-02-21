// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/ConicalVolumeBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/ObjTestWriter.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace tt = boost::test_tools;

namespace Acts {

using namespace UnitLiterals;

// Create a test context
GeometryContext tgContext = GeometryContext();

namespace Test {
BOOST_AUTO_TEST_SUITE(Volumes)

BOOST_AUTO_TEST_CASE(ConcicalVolumeBoundsTests) {
  std::vector<IdentifiedPolyderon> tPolyhedrons;

  auto combineAndDecompose = [&](const SurfacePtrVector& surfaces,
                                 const std::string& name) -> void {
    std::string writeBase = std::string("ConicalVolumeBounds_") + name;

    Polyhedron phCombined;
    size_t is = 0;
    for (const auto& sf : surfaces) {
      Polyhedron phComponent = sf->polyhedronRepresentation(tgContext, 72);
      phCombined += phComponent;
      tPolyhedrons.push_back(
          {writeBase + std::string("_comp_") + std::to_string(is++), false,
           phCombined});
    }
    tPolyhedrons.push_back({writeBase, false, phCombined});
  };

  // Single solid Cone
  ConicalVolumeBounds solidCone(0., 0., 0.45, 50_mm, 50_mm);

  // Test correct parameter return
  BOOST_TEST(solidCone.innerAlpha() == 0.);
  BOOST_TEST(solidCone.innerOffsetZ() == 0.);
  BOOST_TEST(solidCone.outerAlpha() == 0.45);
  BOOST_TEST(solidCone.outerOffsetZ() == 50.);
  BOOST_TEST(solidCone.halflengthZ() == 50.);
  BOOST_TEST(solidCone.averagePhi() == 0.);
  BOOST_TEST(solidCone.halfPhiSector() == M_PI);
  // Derived quantities
  BOOST_TEST(solidCone.innerTanAlpha() == 0.);
  BOOST_TEST(solidCone.innerRmin() == 0.);
  BOOST_TEST(solidCone.innerRmax() == 0.);
  BOOST_TEST(solidCone.outerTanAlpha() == std::tan(0.45));

  double outerRmax = 100_mm * solidCone.outerTanAlpha();
  BOOST_TEST(solidCone.outerRmin() == 0.);
  BOOST_TEST(solidCone.outerRmax() == outerRmax);

  auto solidConeSurfaces = solidCone.decomposeToSurfaces();
  BOOST_TEST(solidConeSurfaces.size() == 2);
  combineAndDecompose(solidConeSurfaces, "Solid");

  // Single solid Cone
  ConicalVolumeBounds cutOffCone(0., 0., 0.45, 80_mm, 50_mm);
  auto cutOffConeSurfaces = cutOffCone.decomposeToSurfaces();
  BOOST_TEST(cutOffConeSurfaces.size() == 3);
  combineAndDecompose(cutOffConeSurfaces, "Cutoff");

  ObjTestWriter::writeObj(tPolyhedrons);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
