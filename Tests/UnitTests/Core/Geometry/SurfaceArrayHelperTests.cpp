// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <random>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/SurfaceArrayHelper.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

using namespace UnitLiterals;

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

using SurfacePtr = std::shared_ptr<const Surface>;
using SurfacePtrVector = std::vector<SurfacePtr>;

// This is what we are testing
SurfaceArrayHelper::Config sahConfig;
SurfaceArrayHelper sArrayHelper = SurfaceArrayHelper(sahConfig);

// Create the surfaces for cylinders and discs
static SurfacePtrVector createCylinder(double r, double halfZ,
                                       unsigned int nPhi, unsigned int nZ,
                                       double rStagger = 0.) {
  // Local y is global z
  Vector3D colY(0., 0., 1.);

  // Prepare module dimensions
  double deltaPhi = 2 * M_PI / nPhi;
  double modHalfY = halfZ / nZ;
  double modHalfX = r * std::tan(deltaPhi * 0.5);
  auto mBounds = std::make_shared<RectangleBounds>(modHalfX, modHalfY);

  SurfacePtrVector cylinderSurfaces;
  for (unsigned int iz = 0; iz < nZ; ++iz) {
    double mz = -halfZ + (2 * iz + 1) * modHalfY;
    double mr = r + (iz % 2) * rStagger;
    for (unsigned int iphi = 0; iphi < nPhi; ++iphi) {
      double mphi = -M_PI + (0.5 * iphi) * deltaPhi;
      Vector3D mCenter(mr * std::cos(mphi), mr * std::sin(mphi), mz);
      Vector3D colZ = mCenter.normalized();
      Vector3D colX = colY.cross(colZ);
      // Create the RotationMatrix
      RotationMatrix3D mRotation;
      mRotation.col(0) = colX;
      mRotation.col(1) = colY;
      mRotation.col(2) = colZ;
      // Get the moduleTransform
      std::shared_ptr<Transform3D> mTransform =
          std::make_shared<Transform3D>(Translation3D(mCenter) * mRotation);
      // Create the surface
      cylinderSurfaces.push_back(
          Surface::makeShared<const PlaneSurface>(mTransform, mBounds));
    }
  }
  return cylinderSurfaces;
}

/// @brief Unit test to see if surfaces are correctly
/// sorted into a cylinder
///
BOOST_AUTO_TEST_CASE(SingleCylinderArray) {
  double r = 33_mm;
  double hZ = 100_mm;
  double stagger = 2_mm;
  unsigned int nPhi = 18, nZ = 12;

  // Create the surfaces
  auto surfaces = createCylinder(r, hZ, nPhi, nZ, stagger);
  BOOST_CHECK(surfaces.size() == nPhi * nZ);

  // Options to sort into one cylinder
  SurfaceArrayHelper::Options sOptions;
  sOptions.rSplitTolerance = 3_mm;

  auto sSurfaces = sArrayHelper.cylinders(
      tgContext, unpack_shared_vector(surfaces), sOptions);

  // This should lead to only one cylinder
  BOOST_CHECK(sSurfaces.size() == 1);
  // All surfaces should be contained
  BOOST_CHECK(sSurfaces[0].surfaces.size() == nPhi * nZ);
  // Check if the reference radius is within tolerance
  BOOST_CHECK(std::abs(sSurfaces[0].splitValue - r) < sOptions.rSplitTolerance);
}

/// @brief Unit test to see if surfaces are correctly
/// sorted into three cylinders
///
BOOST_AUTO_TEST_CASE(ThreeCylinderArray) {
  double r0 = 33_mm, hz0 = 100_mm;
  unsigned int nPhi0 = 18, nZ0 = 12;
  double stagger = 2_mm;

  double r1 = 45_mm, hz1 = 100_mm;
  unsigned int nPhi1 = 38, nZ1 = 12;

  double r2 = 76_mm, hz2 = 100_mm;
  unsigned int nPhi2 = 52, nZ2 = 12;

  // Create the surfaces for the cylinders
  auto surfaces0 = createCylinder(r0, hz0, nPhi0, nZ0, stagger);
  auto surfaces1 = createCylinder(r1, hz1, nPhi1, nZ1, stagger);
  auto surfaces2 = createCylinder(r2, hz2, nPhi2, nZ2, stagger);

  auto surfaces = surfaces0;
  surfaces.insert(surfaces.end(), surfaces1.begin(), surfaces1.end());
  surfaces.insert(surfaces.end(), surfaces2.begin(), surfaces2.end());

  // Check if all surfaces are filled in
  BOOST_CHECK(surfaces.size() == (nPhi0 * nZ0 + nPhi1 * nZ1 + nPhi2 * nZ2));

  // Let's shuffle the surfaces
  std::random_device rd;
  std::mt19937 g(rd());

  std::shuffle(surfaces.begin(), surfaces.end(), g);

  // Options to sort into one cylinder
  SurfaceArrayHelper::Options sOptions;
  sOptions.rSplitTolerance = 3_mm;

  auto sSurfaces = sArrayHelper.cylinders(
      tgContext, unpack_shared_vector(surfaces), sOptions);

  // This should lead into three cylinders
  BOOST_CHECK(sSurfaces.size() == 3);
  // All surfaces should be contained in the right patch
  BOOST_CHECK(sSurfaces[0].surfaces.size() == nPhi0 * nZ0);
  BOOST_CHECK(sSurfaces[1].surfaces.size() == nPhi1 * nZ1);
  BOOST_CHECK(sSurfaces[2].surfaces.size() == nPhi2 * nZ2);
}

}  // namespace Test

}  // namespace Acts