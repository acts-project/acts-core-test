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
          {writeBase + std::string("_comp_") + std::to_string(is++), true,
           phComponent});
    }
    tPolyhedrons.push_back({writeBase, true, phCombined});
  };

  // Single solid Cone
  ConicalVolumeBounds solidCone(0., 0., 0.45, 50_mm, 50_mm, 0., M_PI);

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

  // Single solid Cone - with cut off
  ConicalVolumeBounds cutOffCone(0., 0., 0.45, 80_mm, 50_mm, 0., M_PI);
  auto cutOffConeSurfaces = cutOffCone.decomposeToSurfaces();
  BOOST_TEST(cutOffConeSurfaces.size() == 3);
  combineAndDecompose(cutOffConeSurfaces, "Cutoff");

  // Cone - conical inlay
  ConicalVolumeBounds cutOffHollowCone(0.35, 70_mm, 0.45, 80_mm, 50_mm, 0.,
                                       M_PI);
  auto cutOffHollowConeSurfaces = cutOffHollowCone.decomposeToSurfaces();
  BOOST_TEST(cutOffHollowConeSurfaces.size() == 4);
  combineAndDecompose(cutOffHollowConeSurfaces, "ConicalHollow");

  // Sectoral Cone - conical inlay
  ConicalVolumeBounds cutOffHollowSectoralCone(0.35, 70_mm, 0.45, 80_mm, 50_mm, 0.,
                                       0.456);
  auto cutOffHollowSectoralConeSurfaces = cutOffHollowSectoralCone.decomposeToSurfaces();
  BOOST_TEST(cutOffHollowSectoralConeSurfaces.size() == 6);
  combineAndDecompose(cutOffHollowSectoralConeSurfaces, "SectoralConicalHollow");
  // This one is difficult - let's check inside/outside here
  std::vector<Vector3D> insidePoints; insidePoints.reserve(5000);
  std::vector<Vector3D> outsidePoints; outsidePoints.reserve(5000);
  for (unsigned int iz = 0; iz < 50; ++iz){
    double z = -70_mm + iz*140./50.;
    for (unsigned int ir = 0; ir < 20; ++ir){
      double r = 5_mm + ir * 3_mm;
      for (unsigned int iph = 0; iph < 20; ++iph){
        double phi = -0.6 + iph * 0.06;
        Vector3D pos(r*std::cos(phi),r*std::sin(phi),z);
        if (cutOffHollowSectoralCone.inside(pos,0.)){
          insidePoints.push_back(pos);
        } else {
          outsidePoints.push_back(pos);
        }
      }
    }
  }
  ObjTestWriter::writeObj("ConicalVolumeBounds_inside", insidePoints);
  ObjTestWriter::writeObj("ConicalVolumeBounds_outside", outsidePoints);

  // Single Hollow Cone - cylindrical inlay
  ConicalVolumeBounds cutOffHollowCylCone(10_mm, 0.45, 80_mm, 50_mm, 0., M_PI);
  auto cutOffHollowCylConeSurfaces = cutOffHollowCylCone.decomposeToSurfaces();
  BOOST_TEST(cutOffHollowCylConeSurfaces.size() == 4);
  combineAndDecompose(cutOffHollowCylConeSurfaces, "ConeInnerCylinderHollow");

  // Single Hollow Cylinder - conical inlay
  ConicalVolumeBounds cutOffHollowConeCyl(120_mm, 0.35, 70_mm, 50_mm, 0., M_PI);
  auto cutOffHollowConeCylSurfaces = cutOffHollowConeCyl.decomposeToSurfaces();
  BOOST_TEST(cutOffHollowConeCylSurfaces.size() == 4);
  combineAndDecompose(cutOffHollowConeCylSurfaces, "CylinderInnerConeHollow");
  ObjTestWriter::writeObj(tPolyhedrons);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
