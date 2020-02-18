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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/Units.hpp"

#include <fstream>

namespace Acts {
namespace Test {
namespace Layers {

// Create a test context
GeometryContext tgContext = GeometryContext();

using namespace UnitLiterals;

using SurfaceVector = std::vector<std::shared_ptr<const Surface>>;

SurfaceArrayCreator saCreator;

static void objOutput(const SurfaceVector& sVector, const std::string& name) {
  std::ofstream ostream;
  ostream.open(name + ".obj");
  ObjHelper objH;
  for (const auto& sf : sVector) {
    sf->polyhedronRepresentation(tgContext, 1).draw(objH);
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

BOOST_AUTO_TEST_CASE(CompoundLayerTests) {
  auto rBounds = std::make_shared<RectangleBounds>(10_mm, 4_mm);
  auto bSurfaces =
      barrelSurfaces(32_mm, 1_mm, 26, {-20.5_mm, 0_mm, 20.5_mm}, 0_mm, rBounds);
  objOutput(bSurfaces, "BarrelSurfaces");

  auto tBounds = std::make_shared<TrapezoidBounds>(3_mm, 5_mm, 10_mm);
  auto cRingSurfaces = coneRingSurfaces(26_mm, 1_mm, 26, -40_mm, 0.94, tBounds);
  objOutput(cRingSurfaces, "ConeSurfaces");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Layers
}  // namespace Test
}  // namespace Acts