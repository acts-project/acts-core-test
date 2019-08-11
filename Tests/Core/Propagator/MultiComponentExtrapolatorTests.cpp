// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE MultiComponent Extrapolator Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/MultiMaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/MultiEigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Global definitions
// The path limit abort
using path_limit = detail::PathLimitReached;

std::vector<std::unique_ptr<const Surface>> stepState;

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// Get the navigator and provide the TrackingGeometry
Navigator navigator(tGeometry);
Navigator multi_navigator(tGeometry);

using BFieldType = ConstantBField;
using EigenStepperType = EigenStepper<BFieldType>;
using EigenPropagatorType = Propagator<EigenStepperType, Navigator>;
using MultiEigenStepperType = MultiEigenStepper<BFieldType>;
using MultiEigenPropagatorType = Propagator<MultiEigenStepperType, Navigator>;
using Covariance = BoundSymMatrix;

const double Bz = 2. * units::_T;
BFieldType bField(0, 0, Bz);

// define the scs
EigenStepperType estepper(bField);
EigenPropagatorType epropagator(std::move(estepper), std::move(navigator));

// define the mcs
MultiEigenStepperType multi_estepper(bField);
MultiEigenPropagatorType multi_epropagator(std::move(multi_estepper),
                                           std::move(multi_navigator));

const int ntests = 30;
bool debugMode = false;

// A plane selector for the SurfaceCollector
struct PlaneSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Surface& sf) const {
    return (sf.type() == Surface::Plane);
  }
};

// A direction changing actor
// in MultiMaterialInteractor, at each surface the component splits into 2
// equal ones,
// in this actor, the first component will flip the direction, and it will be
// killed in the MultiStepper since it is a dead component
struct DirChangeActor {
  struct this_result {
    std::vector<Vector3D> dirVec;
  };
  using result_type = this_result;

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    if (state.navigation.currentSurface) {
      if (state.stepping.stateCol.size() > 1) {
        std::get<0>(state.stepping.stateCol.front()).dir =
            -1 * stepper.direction(state.stepping);
        std::get<1>(state.stepping.stateCol.front()) *= 0.0001;
        result.dirVec.push_back(
            std::get<0>(state.stepping.stateCol.front()).dir);
      }
      stepper.normalizeComponents(state.stepping);
    }
  }
};

// This test case checks that no segmentation fault appears in the mcs
// the basic multi stepper propagate
BOOST_DATA_TEST_CASE(
    test_mcs_extrapolation_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution = std::uniform_real_distribution<>(
                       0.4 * units::_GeV, 10. * units::_GeV))) ^
        bdata::random((bdata::seed = 1,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 2,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 3,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 4,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  auto covPtr = std::make_unique<const Covariance>(cov);
  CurvilinearParameters start(std::move(covPtr), pos, mom, q, time);

  using DebugOutput = detail::DebugOutputActor;

  PropagatorOptions<ActionList<DebugOutput>> options(tgContext, mfContext);
  options.debug = debugMode;
  options.maxStepSize = 10. * units::_cm;
  options.pathLimit = 25 * units::_cm;

  const auto& result = multi_epropagator.propagate(start, options).value();
  if (debugMode) {
    const auto& output = result.get<DebugOutput::result_type>();
    std::cout << ">>> mcs Extrapolation output " << std::endl;
    std::cout << output.debugString << std::endl;
  }
  BOOST_CHECK(result.endParameters != nullptr);
}

// This test case checks that no segmentation fault appears
// - this tests the same surfaceHit of different stepper
// - this tests creat 4 propagation stream
// 1. single stepper
// 2. multi stepper with 1 component
// 3. multi stepper with multi-material-interactor : when meet materail
// surface, each component copy themselves into 2 components
// 4. multi stepper with MultiMaterialInteractor and a DiractionChange actor
// : when meet materail surface, each component copy itselves into 2
// components, then the first component flip its direction and reweight it to
// a small number to mimic a dead component
// all the above cases should collect same surface hit results

BOOST_DATA_TEST_CASE(
    test_equal_scs_mcs_collection_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution = std::uniform_real_distribution<>(
                       0.4 * units::_GeV, 10. * units::_GeV))) ^
        bdata::random((bdata::seed = 1,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 2,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 3,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 14,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  auto covPtr = std::make_unique<const Covariance>(cov);
  CurvilinearParameters start(std::move(covPtr), pos, mom, q, time);

  // A PlaneSelector for the SurfaceCollector
  using PlaneCollector = SurfaceCollector<PlaneSelector>;

  PropagatorOptions<ActionList<PlaneCollector>> options(tgContext, mfContext);
  options.maxStepSize = 10. * units::_cm;
  options.pathLimit = 25 * units::_cm;
  options.debug = debugMode;

  PropagatorOptions<ActionList<PlaneCollector>> multi_options(tgContext,
                                                              mfContext);
  multi_options.maxStepSize = 10. * units::_cm;
  multi_options.pathLimit = 25 * units::_cm;
  multi_options.debug = debugMode;

  PropagatorOptions<ActionList<PlaneCollector, MultiMaterialInteractor>>
      multi_material_options(tgContext, mfContext);
  multi_material_options.maxStepSize = 10. * units::_cm;
  multi_material_options.pathLimit = 25 * units::_cm;
  multi_material_options.debug = debugMode;

  PropagatorOptions<
      ActionList<PlaneCollector, MultiMaterialInteractor, DirChangeActor>>
      flip_options(tgContext, mfContext);
  flip_options.maxStepSize = 10. * units::_cm;
  flip_options.pathLimit = 25 * units::_cm;
  flip_options.debug = debugMode;

  // sigle component
  const auto& result = epropagator.propagate(start, options).value();
  auto collector_result = result.get<PlaneCollector::result_type>();

  // multi component
  const auto& multi_result =
      multi_epropagator.propagate(start, multi_options).value();
  auto multi_collector_result = multi_result.get<PlaneCollector::result_type>();

  // multi_material component
  const auto& multi_material_result =
      multi_epropagator.propagate(start, multi_material_options).value();
  auto multi_material_collector_result =
      multi_material_result.get<PlaneCollector::result_type>();

  // flip the multi component
  const auto& flip_result =
      multi_epropagator.propagate(start, flip_options).value();
  auto flip_collector_result = flip_result.get<PlaneCollector::result_type>();

  // check if the size of collected results are equal
  BOOST_CHECK_EQUAL(collector_result.collected.size(),
                    multi_collector_result.collected.size());
  BOOST_CHECK_EQUAL(collector_result.collected.size(),
                    multi_material_collector_result.collected.size());
  BOOST_CHECK_EQUAL(collector_result.collected.size(),
                    flip_collector_result.collected.size());

  // check if multi-stepper with 1 component collect the same result with
  // single-stepper
  BOOST_CHECK(collector_result.collected == multi_collector_result.collected);
  // check if multi-stepper with several (by material split) components
  // collect the same result with single-stepper
  BOOST_CHECK(collector_result.collected ==
              multi_material_collector_result.collected);
  // check if multi-stepper with several (by material split, but each time
  // kill the flip one) components collect the same result with single-stepper
  BOOST_CHECK(collector_result.collected == flip_collector_result.collected);
}

}  // namespace Test
}  // namespace Acts
