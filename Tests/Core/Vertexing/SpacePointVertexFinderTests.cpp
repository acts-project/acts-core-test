// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE SpacePointVertexFinder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <climits>
#include <memory>
#include <random>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Vertexing/SpacePointVertexFinder.hpp"
#include "Acts/Vertexing/SpacePointVertexScanner.hpp"
#include "Acts/Vertexing/SpacePointVertexProjectors.hpp"
#include "Acts/Vertexing/SpacePointVertexTargets.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Global definitions
using EndOfWorld = detail::EndOfWorldReached;

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// create a navigator for this tracking geometry
Navigator navigator(tGeometry);

using BField = ConstantBField;
using EigenStepper = EigenStepper<BField>;
using EigenPropagator = Propagator<EigenStepper, Navigator>;
using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

StraightLineStepper slstepper;
StraightLinePropagator slpropagator(std::move(slstepper), std::move(navigator));

std::vector<double> tvertices = {-5_mm, 2_mm, 0.3_mm, 10_mm, 25_mm};
const unsigned int nTracks = 100;
const double thetaH = 0.35;

std::uniform_real_distribution<> uniform(0., 1.);
std::default_random_engine generator(42);

BOOST_AUTO_TEST_CASE(cylindrical_target_test) {
  struct TestPoint {
    Eigen::Array<double, 1, 3> rawValues;
  };

  // Create a cylinder target
  Eigen::Array<double, 2, 2> minmaxc;
  minmaxc << -400., 400., 30, 36.;

  CylindricalDetectorTarget cTarget(minmaxc, 5);
  BOOST_TEST(cTarget.targetType() == 0);
  BOOST_TEST(cTarget.targetID() == 5);
  BOOST_TEST(double(cTarget) == 33.);

  // Get a point on the cylinder
  Vector3D centrallyOn(33., 0., 0.);
  TestPoint cop;
  cop.rawValues = cTarget.rawValues(centrallyOn);
  BOOST_TEST(cop.rawValues(0, 2) == 33.);

  Vector3D midneg(34., 0., -200.);
  TestPoint min;
  min.rawValues = cTarget.rawValues(midneg);
  auto ppars = cTarget.projectionParameters(min);
  BOOST_TEST(ppars(0, 1) == -200.);
  cTarget.scaleProjection(ppars);
  BOOST_TEST(ppars(0, 1) == -0.5);

  // Create a disk target - at negative side
  Eigen::Array<double, 2, 2> minmaxd;
  minmaxd << -415., -405., 40, 80.;
  CylindricalDetectorTarget dTargetNeg(minmaxd, 6);
  BOOST_TEST(dTargetNeg.targetType() == 1);
  BOOST_TEST(dTargetNeg.targetID() == 6);
  BOOST_TEST(double(dTargetNeg) == -410.);
  Vector3D midr(60, 0., -408.);
  TestPoint mir;
  mir.rawValues = dTargetNeg.rawValues(midr);
  ppars = dTargetNeg.projectionParameters(mir);
  BOOST_TEST(ppars(0, 1) == 60.);
  dTargetNeg.scaleProjection(ppars);
  BOOST_TEST(ppars(0, 1) == -0.5);

  // Create a disk target - at positive side
  Eigen::Array<double, 2, 2> minmaxp;
  minmaxd << 405., 415., 40, 80.;
  CylindricalDetectorTarget dTargetPos(minmaxd, 7);
  BOOST_TEST(dTargetPos.targetType() == 1);
  BOOST_TEST(dTargetPos.targetID() == 7);
  BOOST_TEST(double(dTargetPos) == 410.);
  midr = Vector3D(60, 0., 408.);
  mir.rawValues = dTargetPos.rawValues(midr);
  ppars = dTargetPos.projectionParameters(mir);
  BOOST_TEST(ppars(0, 1) == 60.);
  dTargetPos.scaleProjection(ppars);
  BOOST_TEST(ppars(0, 1) == 0.5);
}

BOOST_AUTO_TEST_CASE(straight_tracks_test) {
  // Action list for direct navigator with its initalizer
  using ActionList = ActionList<SurfaceCollector<>>;
  using AbortList = AbortList<EndOfWorld>;

  // Direct options definition
  using Options = PropagatorOptions<ActionList, AbortList>;
  Options pOptions(tgContext, mfContext);
  // Surface collector configuration
  auto& hCollector = pOptions.actionList.template get<SurfaceCollector<>>();
  hCollector.selector.selectSensitive = true;

  struct Hit {
    Vector3D xyz;
    unsigned int lid = std::numeric_limits<unsigned int>::max();

    // Return the poistion
    const Vector3D& position() const { return xyz; }

    // Return the layer id
    unsigned int layer() const { return lid; }
  };

  // Create the space points per layer
  using HitVector = std::vector<Hit>;
  std::array<HitVector, 3> seeds;
  HitVector confirmations;

  using ConfirmationTargets = std::vector<CylindricalDetectorTarget>;

  using Scanner =
      SpacePointVertexScanner<360, LinearProjector, ConfirmationTargets>;
  Scanner::Config scannerCfg;
  scannerCfg.confirmBins0 = 720;
  scannerCfg.confirmBins1 = 1000;
  scannerCfg.neighborMask = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1},
                             {0, 1},   {1, -1}, {1, 0},  {1, 1}};
  Scanner scanner(scannerCfg);

  using VertexFinder = SpacePointVertexFinder<Scanner>;
  VertexFinder::Config finderConfig;
  finderConfig.targetBins = 400;
  finderConfig.targetRange = {-25., 25.};
  VertexFinder finder(finderConfig, std::move(scanner));

  // Loop over every vertex
  for (const auto& tvz : tvertices) {
    // Create nTracks track per vertex
    for (unsigned int itrack = 0; itrack < nTracks; ++itrack) {
      double pT = 100_GeV;
      double phi = M_PI * (2 * uniform(generator) - 1);
      double theta = 0.5 * M_PI + (thetaH * (2 * uniform(generator) - 1));
      double charge = uniform(generator) > 0.5 ? 1. : -1.;

      // define start parameters
      double x = 0;
      double y = 0;
      double z = tvz;
      double px = pT * cos(phi);
      double py = pT * sin(phi);
      double pz = pT / tan(theta);
      double q = charge;
      double stime = 0.;
      Vector3D pos(x, y, z);
      Vector3D mom(px, py, pz);
      CurvilinearParameters start(std::nullopt, pos, mom, q, stime);

      auto pResult = slpropagator.propagate(start, pOptions);
      if (pResult.ok()) {
        const auto& pOutput = pResult.value();
        auto& cSurfaces =
            pOutput.template get<SurfaceCollector<>::result_type>();
        // Collect the hits
        for (const auto& cSurface : cSurfaces.collected) {
          // Get the layer ID
          unsigned int layerID = cSurface.surface->geoID().layer() / 2;
          if (layerID < 4) {
            seeds[layerID - 1].push_back(Hit{cSurface.position, layerID});
          } else {
            confirmations.push_back(Hit{cSurface.position, layerID});
          }
        }
      }
    }
  }

  using seedComb = std::vector<unsigned int>;
  std::vector<seedComb> seedCombinations = {{0, 1}, {0, 2}, {1, 2}};
  std::vector<HitVector> confirmationContainer = {confirmations};
  auto hist = finder.find(seeds, seedCombinations, confirmationContainer);
}

}  // namespace Test
}  // namespace Acts