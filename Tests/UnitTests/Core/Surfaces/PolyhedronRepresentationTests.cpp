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
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"

// Cone surface
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"

// Cylinder surface
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

// Disc Surface
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

// Plane Surface
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/Units.hpp"

#include <fstream>

namespace Acts {

using namespace UnitLiterals;

using IdentifiedPolyderon = std::pair<std::string, PolyhedronRepresentation>;

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

/// Helper method to write
static void writeSectorLinesObj(const std::string& name,
                                const std::pair<Vector3D, Vector3D>& lineA,
                                const std::pair<Vector3D, Vector3D>& lineB) {
  std::ofstream ostream;
  ostream.open(name + ".obj");
  ObjHelper objH;
  objH.line(lineA.first, lineA.second);
  objH.line(lineB.first, lineB.second);
  objH.write(ostream);
  ostream.close();
}

/// Helper method to write
static void writeSectorPlanesObj(const std::string& name, double phiSec,
                                 double phiAvg, double hX, double hY) {
  // Construct the helper planes for sectoral building
  auto sectorBounds = std::make_shared<RectangleBounds>(hX, hY);

  Vector3D helperColX(0., 0., 1.);
  Vector3D helperColY(1., 0., 0.);
  Vector3D helperColZ(0., 1., 0.);
  RotationMatrix3D helperRotation;
  helperRotation.col(0) = helperColX;
  helperRotation.col(1) = helperColY;
  helperRotation.col(2) = helperColZ;
  // curvilinear surfaces are boundless
  Transform3D helperTransform{helperRotation};

  auto sectorTransformM = std::make_shared<Transform3D>(helperTransform);
  sectorTransformM->prerotate(AngleAxis3D(phiAvg - phiSec, helperColX));

  auto sectorTransformP = std::make_shared<Transform3D>(helperTransform);
  sectorTransformP->prerotate(AngleAxis3D(phiAvg + phiSec, helperColX));

  auto sectorPlaneM =
      Surface::makeShared<PlaneSurface>(sectorTransformM, sectorBounds);

  auto sectorPlaneP =
      Surface::makeShared<PlaneSurface>(sectorTransformP, sectorBounds);

  std::ofstream ostream;
  ostream.open(name + ".obj");
  ObjHelper objH;
  sectorPlaneM->polyhedronRepresentation(tgContext).draw(objH);
  sectorPlaneP->polyhedronRepresentation(tgContext).draw(objH);
  objH.write(ostream);
  ostream.close();
}

/// Helper method to be called from sub tests
static void writeObj(const std::vector<IdentifiedPolyderon>& iphs) {
  for (const auto& iph : iphs) {
    std::ofstream ostream;
    ostream.open(iph.first + ".obj");
    ObjHelper objH;
    iph.second.draw(objH);
    objH.write(ostream);
    ostream.close();
  }
}

std::vector<std::pair<std::string, unsigned int>> testModes = {
    {"", 72}, {"Extremas", 1}, {"TriangleMesh", 72}};

auto transform = std::make_shared<Transform3D>(Transform3D::Identity());

BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(PolyhedronRepresentationTest) {
  std::vector<Vector3D> vertices = {Vector3D(-1, -1, 0.), Vector3D(1., -1, 0.),
                                    Vector3D(0., 1., 0.)};
  std::vector<std::vector<size_t>> faces = {{0, 1, 2}};

  PolyhedronRepresentation triangle(vertices, faces);
  BOOST_CHECK(vertices == triangle.vertices);
  BOOST_CHECK(faces == triangle.faces);
}

/// Unit tests for Cone Surfaces
BOOST_AUTO_TEST_CASE(ConeSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  double hzpos = 35_mm;
  double hzneg = -20_mm;
  double alpha = 0.234;
  double phiSector = 0.358;
  writeSectorPlanesObj("ConeSectorPlanes", phiSector, 0., hzpos, hzpos);

  for (const auto& mode : testModes) {
    // For ConeSurface types the polyhedron representation is
    /// always a triangular mesh
    bool triangulate = (mode.first == "TriangleMesh");
    if (triangulate) {
      continue;
    }

    /// The full cone on one side
    auto cone = std::make_shared<ConeBounds>(alpha, 0_mm, hzpos);
    auto oneCone = Surface::makeShared<ConeSurface>(transform, cone);
    auto oneConePh = oneCone->polyhedronRepresentation(tgContext, mode.second);
    size_t expectedFaces = mode.second < 4 ? 4 : mode.second;
    BOOST_CHECK_EQUAL(oneConePh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(oneConePh.vertices.size(), expectedFaces + 1);
    // Check the extent in space
    double r = hzpos * std::tan(alpha);
    auto extent = oneConePh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.xrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hzpos, 1e-6);
    testTypes.push_back({"ConeOneFull" + mode.first, oneConePh});

    // The full cone on both sides
    auto coneBoth = std::make_shared<ConeBounds>(alpha, hzneg, hzpos);
    auto twoCones = Surface::makeShared<ConeSurface>(transform, coneBoth);
    auto twoConesPh =
        twoCones->polyhedronRepresentation(tgContext, mode.second);
    expectedFaces = mode.second < 4 ? 8 : 2 * mode.second;
    BOOST_CHECK_EQUAL(twoConesPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(twoConesPh.vertices.size(), expectedFaces + 1);
    extent = twoConesPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.xrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, hzneg, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hzpos, 1e-6);
    testTypes.push_back({"ConesTwoFull" + mode.first, twoConesPh});

    // A centered sectoral cone on both sides
    auto sectoralBoth =
        std::make_shared<ConeBounds>(alpha, hzneg, hzpos, phiSector, 0.);
    auto sectoralCones =
        Surface::makeShared<ConeSurface>(transform, sectoralBoth);
    auto sectoralConesPh =
        sectoralCones->polyhedronRepresentation(tgContext, mode.second);
    extent = sectoralConesPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, hzneg, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hzpos, 1e-6);
    testTypes.push_back({"ConesSectoral" + mode.first, sectoralConesPh});
  }
  writeObj(testTypes);
}

/// Unit tests for Cylinder Surfaces
BOOST_AUTO_TEST_CASE(CylinderSurfacePolyhedrons) {
  double r = 25_mm;
  double hZ = 35_mm;

  double phiSector = 0.458;
  double averagePhi = -1.345;
  writeSectorPlanesObj("CylinderCentralSectorPlanes", phiSector, 0., 1.5 * r,
                       1.5 * hZ);
  writeSectorPlanesObj("CylinderShiftedSectorPlanes", phiSector, averagePhi,
                       1.5 * r, 1.5 * hZ);

  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    bool triangulate = (mode.first == "TriangleMesh");
    size_t mulitplier = triangulate ? 2 : 1;

    size_t expectedFaces = mode.second < 4 ? 4 : mulitplier * mode.second;
    size_t expectedVertices = mode.second < 4 ? 8 : 2 * mode.second;

    /// The full cone on one side
    auto cylinder = std::make_shared<CylinderBounds>(r, hZ);
    auto fullCylinder =
        Surface::makeShared<CylinderSurface>(transform, cylinder);
    auto fullCylinderPh = fullCylinder->polyhedronRepresentation(
        tgContext, mode.second, triangulate);

    BOOST_CHECK_EQUAL(fullCylinderPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(fullCylinderPh.vertices.size(), expectedVertices);
    // Check the extent in space
    auto extent = fullCylinderPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.xrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hZ, 1e-6);
    testTypes.push_back({"CylinderFull" + mode.first, fullCylinderPh});

    /// The full cone on one side
    auto sectorCentered = std::make_shared<CylinderBounds>(r, phiSector, hZ);
    auto centerSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorCentered);
    auto centerSectoredCylinderPh =
        centerSectoredCylinder->polyhedronRepresentation(tgContext, mode.second,
                                                         triangulate);

    // Check the extent in space
    extent = centerSectoredCylinderPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.first, r * std::cos(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.xrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.first, -r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.second, r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hZ, 1e-6);
    testTypes.push_back(
        {"CylinderSectorCentered" + mode.first, centerSectoredCylinderPh});

    /// The full cone on one side
    auto sectorShifted =
        std::make_shared<CylinderBounds>(r, averagePhi, phiSector, hZ);
    auto shiftedSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorShifted);
    auto shiftedSectoredCylinderPh =
        shiftedSectoredCylinder->polyhedronRepresentation(
            tgContext, mode.second, triangulate);

    // Check the extent in space
    extent = shiftedSectoredCylinderPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.rrange.first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, hZ, 1e-6);
    testTypes.push_back(
        {"CylinderSectorShifted" + mode.first, shiftedSectoredCylinderPh});
  }

  writeObj(testTypes);
}

/// Unit tests for Disc Surfaces
BOOST_AUTO_TEST_CASE(DiscSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  double innerR = 10_mm;
  double outerR = 25_mm;
  double phiSector = 0.345;
  double averagePhi = -2.0;

  double cphi = std::cos(phiSector);
  double sphi = std::sin(phiSector);

  std::pair<Vector3D, Vector3D> lineA = {
      Vector3D(0., 0., 0.), Vector3D(outerR * cphi, outerR * sphi, 0.)};
  std::pair<Vector3D, Vector3D> lineB = {
      Vector3D(0., 0., 0.), Vector3D(outerR * cphi, -outerR * sphi, 0.)};
  writeSectorLinesObj("DiscSectorLines", lineA, lineB);

  double minPhi = averagePhi - phiSector;
  double maxPhi = averagePhi + phiSector;
  linaA = {Vector3D(0., 0., 0.),
           Vector3D(outerR * std::cos(minPhi), outerR * std::sin(minPhi), 0.)};
  linaB = {Vector3D(0., 0., 0.),
           Vector3D(outerR * std::cos(maxPhi), outerR * std::sin(maxPhi), 0.)};
  writeSectorLinesObj("DiscSectorLinesShifted", lineA, lineB);

  for (const auto& mode : testModes) {
    bool triangulate = (mode.first == "TriangleMesh");

    // Full disc
    auto disc = std::make_shared<RadialBounds>(0_mm, outerR);
    auto fullDisc = Surface::makeShared<DiscSurface>(transform, disc);
    auto fullDiscPh =
        fullDisc->polyhedronRepresentation(tgContext, mode.second, triangulate);

    unsigned int expectedVertices = mode.second > 4 ? mode.second : 4;
    unsigned int expectedFaces = triangulate ? 72 : 1;
    if (triangulate) {
      expectedVertices += 1;
    }
    BOOST_CHECK_EQUAL(fullDiscPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(fullDiscPh.vertices.size(), expectedVertices);

    auto extent = fullDiscPh.surfaceExtent();
    CHECK_CLOSE_ABS(extent.xrange.first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.xrange.second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.yrange.second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.rrange.second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.zrange.second, 0., 1e-6);

    testTypes.push_back({"DiscFull" + mode.first, fullDiscPh});

    // Ring disc
    auto radial = std::make_shared<RadialBounds>(innerR, outerR);
    auto radialDisc = Surface::makeShared<DiscSurface>(transform, radial);
    auto radialPh = radialDisc->polyhedronRepresentation(tgContext, mode.second,
                                                         triangulate);
    testTypes.push_back({"DiscRing" + mode.first, radialPh});

    // Sectoral disc - around 0.
    auto sector = std::make_shared<RadialBounds>(innerR, outerR, phiSector);
    auto sectorDisc = Surface::makeShared<DiscSurface>(transform, sector);
    auto sectorPh = sectorDisc->polyhedronRepresentation(tgContext, mode.second,
                                                         triangulate);
    testTypes.push_back({"DiscSectorCentered" + mode.first, sectorPh});

    // Sectoral disc - shifted
    auto sectorShifted =
        std::make_shared<RadialBounds>(innerR, outerR, phiSector, averagePhi);
    auto sectorDiscShifted =
        Surface::makeShared<DiscSurface>(transform, sectorShifted);
    auto sectorPhShifted = sectorDiscShifted->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    testTypes.push_back({"DiscSectorShifted" + mode.first, sectorPhShifted});
  }
  writeObj(testTypes);
}

/// Unit tests for Plane Surfaces
BOOST_AUTO_TEST_CASE(PlaneSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    bool triangulate = (mode.first == "TriangleMesh");

    /// Rectangular Plane
    auto rectangular = std::make_shared<RectangleBounds>(10_mm, 25_mm);
    auto rectangularPlane =
        Surface::makeShared<PlaneSurface>(transform, rectangular);
    auto rectangularPh = rectangularPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(rectangularPh.vertices.size() == 4);
    if (not triangulate) {
      BOOST_CHECK(rectangularPh.faces.size() == 1);
      std::vector<size_t> expectedRect = {0, 1, 2, 3};
      BOOST_CHECK(rectangularPh.faces[0] == expectedRect);
    } else {
      BOOST_CHECK(rectangularPh.faces.size() == 2);
    }
    testTypes.push_back({"PlaneRectangle" + mode.first, rectangularPh});

    /// Trapezoidal Plane
    auto trapezoid = std::make_shared<TrapezoidBounds>(10_mm, 25_mm, 35_mm);
    auto trapezoidalPlane =
        Surface::makeShared<PlaneSurface>(transform, trapezoid);
    auto trapezoidalPh = trapezoidalPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(trapezoidalPh.vertices.size() == 4);
    if (not triangulate) {
      BOOST_CHECK(trapezoidalPh.faces.size() == 1);
      std::vector<size_t> expectedTra = {0, 1, 2, 3};
      BOOST_CHECK(trapezoidalPh.faces[0] == expectedTra);
    } else {
      BOOST_CHECK(trapezoidalPh.faces.size() == 2);
    }
    testTypes.push_back({"PlaneTrapezoid" + mode.first, trapezoidalPh});

    /// Full ellispoidal Plane
    /**
    auto ellipse = std::make_shared<EllipseBounds>(10_mm,20_mm,30_mm,40_mm);
    auto ellipsoidPlane = Surface::makeShared<PlaneSurface>(transform,ellipse);
    auto ellispoidPh = ellipsoidPlane->polyhedronRepresentation(tgContext,72);
    BOOST_CHECK(ellispoidPh.vertices.size()==72);
    testTypes.push_back({"PlaneFullEllipse",ellispoidPh});
    */

    /// Sectoral epplise

    /// Diamond shaped plane
    auto diamond =
        std::make_shared<DiamondBounds>(10_mm, 20_mm, 15_mm, 40_mm, 50_mm);
    auto diamondPlane = Surface::makeShared<PlaneSurface>(transform, diamond);
    auto diamondPh = diamondPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(diamondPh.vertices.size() == 6);
    if (not triangulate) {
      BOOST_CHECK(diamondPh.faces.size() == 1);
    } else {
      BOOST_CHECK(diamondPh.faces.size() == 4);
    }
    testTypes.push_back({"PlaneDiamond" + mode.first, diamondPh});
  }
  writeObj(testTypes);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts