// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CompoundLayer.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/Units.hpp"

#include <fstream>
#include <iostream>

namespace Acts {
namespace Test {
namespace Layers {

// Create a test context
GeometryContext tgContext = GeometryContext();

using namespace UnitLiterals;

using SurfaceVector = std::vector<std::shared_ptr<const Surface>>;

SurfaceArrayCreator saCreator;

static void objOutput(const SurfaceVector& sVector, const std::string& name,
                      unsigned int lseg = 1) {
  std::ofstream ostream;
  ostream.open(name + ".obj");
  ObjHelper objH;
  for (const auto& sf : sVector) {
    sf->polyhedronRepresentation(tgContext, lseg).draw(objH);
  }
  objH.write(ostream);
  ostream.close();
}

static SurfaceVector coneRingSurfaces(
    double rref, double rstag, size_t nPhi, double z, double skew,
    std::shared_ptr<const TrapezoidBounds> trdbounds) {
  SurfaceVector res;

  double phiStep = 2 * M_PI / nPhi;
  for (size_t iphi = 0; iphi < nPhi; ++iphi) {
    double r = rref + (iphi % 2) * rstag;
    double phi = std::fma(iphi, phiStep, 0.);
    Transform3D trans = Transform3D::Identity();
    trans.translate(Vector3D(r * std::cos(phi), r * std::sin(phi), z));
    trans.rotate(Eigen::AngleAxisd(phi - 0.5 * M_PI, Vector3D(0, 0, 1)));
    trans.rotate(Eigen::AngleAxisd(skew, Vector3D(1, 0, 0)));
    auto transptr = std::make_shared<const Transform3D>(trans);
    std::shared_ptr<Surface> srf =
        Surface::makeShared<PlaneSurface>(transptr, trdbounds);
    res.push_back(srf);
  }
  return res;
}

static SurfaceVector barrelSurfaces(
    double rref, double rstag, size_t nPhi, std::vector<double> zcenters,
    double zstag, std::shared_ptr<const RectangleBounds> rbounds) {
  SurfaceVector res;

  double phiStep = 2 * M_PI / nPhi;

  for (size_t iz = 0; iz < zcenters.size(); ++iz) {
    double z = zcenters[iz] + (iz % 2) * zstag;
    for (size_t iphi = 0; iphi < nPhi; ++iphi) {
      double r = rref + (iphi % 2) * rstag;
      double phi = std::fma(iphi, phiStep, 0.);
      Transform3D trans = Transform3D::Identity();
      trans.rotate(Eigen::AngleAxisd(phi, Vector3D(0, 0, 1)));
      trans.translate(Vector3D(r, 0, z));
      trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3D(0, 1, 0)));
      auto transptr = std::make_shared<const Transform3D>(trans);
      std::shared_ptr<Surface> srf =
          Surface::makeShared<PlaneSurface>(transptr, rbounds);

      res.push_back(srf);
    }
  }
  return res;
}

BOOST_AUTO_TEST_SUITE(Layers)

/// This test is a step by step test how to build a compound layer
///
BOOST_AUTO_TEST_CASE(CompoundLayerTests) {
  double envelopeR = 2_mm;
  double envelopeZ = 2_mm;

  Extent bExtent;
  auto rBounds = std::make_shared<RectangleBounds>(10_mm, 4_mm);
  auto bSurfaces =
      barrelSurfaces(32_mm, 1_mm, 26, {-20.5_mm, 0_mm, 20.5_mm}, 0_mm, rBounds);
  objOutput(bSurfaces, "BarrelSurfaces");

  for (const auto& bsf : bSurfaces) {
    bExtent += bsf->polyhedronRepresentation(tgContext, 1).extent();
  }

  std::cout << "Dimensions of the Barrel section: " << std::endl;
  std::cout << bExtent << std::endl;

  Extent cExtent;
  std::vector<Vector3D> vertices;
  auto tBounds = std::make_shared<TrapezoidBounds>(3_mm, 5_mm, 10_mm);
  auto cSurfaces = coneRingSurfaces(28_mm, 1_mm, 26, -40_mm, 1.15, tBounds);
  objOutput(cSurfaces, "ConeSurfaces");

  for (const auto& csf : cSurfaces) {
    auto csfph = csf->polyhedronRepresentation(tgContext, 1);
    vertices.insert(vertices.end(), csfph.vertices.begin(),
                    csfph.vertices.end());
    cExtent += csfph.extent();
  }

  std::cout << "Dimensions of the Cone section: " << std::endl;
  std::cout << cExtent << std::endl;

  Extent narrowEnd;
  Extent openEnd;
  for (const auto& v : vertices) {
    if (std::abs(v.z() - cExtent.min(binZ)) < envelopeZ) {
      narrowEnd.check(v);
    } else if (std::abs(v.z() - cExtent.max(binZ)) < envelopeZ) {
      openEnd.check(v);
    }
  }

  std::cout << "Dimensions of the narrow end of the cone: " << std::endl;
  std::cout << narrowEnd << std::endl;

  std::cout << "Dimensions of the open end of the cone: " << std::endl;
  std::cout << openEnd << std::endl;

  // Barrel basic parameters
  double bInnerR = bExtent.min(binR) - envelopeR;
  double bOuterR = bExtent.max(binR) + envelopeR;
  double bHalfZ = 0.5 * bExtent.range(binZ) + envelopeZ;

  // Cone basic parameters
  double cZmin = narrowEnd.min(binZ) - envelopeZ;
  double cZmax = openEnd.max(binZ) + envelopeZ;
  double cDeltaZ = cZmax - cZmin;

  double cInnerRmin = narrowEnd.min(binR) - envelopeR;
  double cInnerRmax = openEnd.min(binR) - envelopeR;
  double cDeltaInnerR = cInnerRmax - cInnerRmin;
  double cInnerAlpha = std::atan2(cDeltaInnerR, cDeltaZ);
  double cInnerZ = cZmin - cInnerRmin / std::tan(cInnerAlpha);

  double cOuterRmin = narrowEnd.max(binR) + envelopeR;
  double cOuterRmax = openEnd.max(binR) + envelopeR;
  double cDeltaOuterR = cOuterRmax - cOuterRmin;
  double cOuterAlpha = std::atan2(cDeltaOuterR, cDeltaZ);
  double cOuterZ = cZmin - cOuterRmin / std::tan(cOuterAlpha);

  // @TODO harmonize
  auto bTransform = std::make_shared<Transform3D>(Transform3D::Identity());

  auto bInnerBounds = std::make_shared<CylinderBounds>(bInnerR, bHalfZ);
  auto bInnerSurface =
      Surface::makeShared<CylinderSurface>(bTransform, bInnerBounds);

  auto bOuterBounds = std::make_shared<CylinderBounds>(bOuterR, bHalfZ);
  auto bOuterSurface =
      Surface::makeShared<CylinderSurface>(bTransform, bOuterBounds);

  auto cInnerTransform = std::make_shared<Transform3D>(Transform3D::Identity());
  cInnerTransform->pretranslate(Vector3D(0., 0., cInnerZ));

  auto cInnerSurface = Surface::makeShared<ConeSurface>(
      cInnerTransform, cInnerAlpha, cZmin - cInnerZ, cZmax - cInnerZ);

  auto cOuterTransform = std::make_shared<Transform3D>(Transform3D::Identity());
  cOuterTransform->pretranslate(Vector3D(0., 0., cOuterZ));

  auto cOuterSurface = Surface::makeShared<ConeSurface>(
      cOuterTransform, cOuterAlpha, cZmin - cOuterZ, cZmax - cOuterZ);

  SurfaceVector bApproachSurfaces = {bInnerSurface, bOuterSurface};
  objOutput(bApproachSurfaces, "BarrelApproachSurfaces", 72);

  SurfaceVector cApproachSurfaces = {cInnerSurface, cOuterSurface};
  objOutput(cApproachSurfaces, "ConeApproachSurfaces", 72);

  SurfaceVector approachSurfaces = {cInnerSurface, bInnerSurface, cOuterSurface,
                                    bOuterSurface};
  objOutput(approachSurfaces, "ApproachSurfaces", 72);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Layers
}  // namespace Test
}  // namespace Acts