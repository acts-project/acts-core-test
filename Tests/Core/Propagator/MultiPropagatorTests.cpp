// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Propagator Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/MultiEigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Utilities/MagneticFieldContext.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

using Acts::VectorHelpers::perp;

namespace Acts {

namespace Test {

  // Create a test context
  GeometryContext      tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();

  using cstep                    = detail::ConstrainedStep;
  using BFieldType               = ConstantBField;
  using EigenStepperType         = EigenStepper<BFieldType>;
  using MultiEigenStepperType    = MultiEigenStepper<BFieldType>;
  using EigenPropagatorType      = Propagator<EigenStepperType>;
  using MultiEigenPropagatorType = Propagator<MultiEigenStepperType>;

  const double          Bz = 2. * units::_T;
  BFieldType            bField(0, 0, Bz);
  EigenStepperType      stepper(bField);
  MultiEigenStepperType mcstepper(bField);

  auto aSurface
      = Surface::makeShared<CylinderSurface>(nullptr, 10., 1000. * units::_mm);
  auto bSurface
      = Surface::makeShared<CylinderSurface>(nullptr, 20., 1000. * units::_mm);
  auto cSurface
      = Surface::makeShared<CylinderSurface>(nullptr, 30., 1000. * units::_mm);

  const int ntests = 100;

  BOOST_DATA_TEST_CASE(
      test_eigenstepper_targetSurface_and_onSurface,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    using ActionListType      = ActionList<>;
    using AbortConditionsType = AbortList<>;
    using PropagateState      = EigenPropagatorType::
        State<PropagatorOptions<ActionListType, AbortConditionsType>>;

    // setup propagation options
    PropagatorOptions<ActionListType, AbortConditionsType> options(tgContext,
                                                                   mfContext);

    // define start parameters
    double                           x  = 0;
    double                           y  = 0;
    double                           z  = 0;
    double                           px = pT * cos(phi);
    double                           py = pT * sin(phi);
    double                           pz = pT / tan(theta);
    double                           q  = dcharge;
    Vector3D                         pos(x, y, z);
    Vector3D                         mom(px, py, pz);
    CurvilinearParameters            start(nullptr, pos, mom, q);
    PropagateState                   prop_state(start, options);
    VoidIntersectionCorrector        corrector;
    Acts::NavigationOptions<Surface> opt(forward, true);

    // target the surface and update the stepSize
    auto intersect1 = stepper.targetSurface(
        prop_state.stepping, aSurface.get(), opt, corrector);
    // calculate with the intersection estimate
    auto intersect2 = aSurface.get()->surfaceIntersectionEstimate(
        prop_state.geoContext,
        stepper.position(prop_state.stepping),
        stepper.direction(prop_state.stepping),
        opt,
        corrector);
    // test if targetSurface find the correct intersection
    CHECK_CLOSE_REL(intersect1.intersection.position,
                    intersect2.intersection.position,
                    10e-6);
    CHECK_CLOSE_REL(intersect1.intersection.pathLength,
                    intersect2.intersection.pathLength,
                    10e-6);

    // propagate to the surface
    stepper.step(prop_state);
    // test if the component is on the surface
    bool isOnSurface
        = stepper.surfaceReached(prop_state.stepping, aSurface.get());
    BOOST_CHECK_EQUAL(isOnSurface, 1);
  }

  BOOST_DATA_TEST_CASE(
      test_multieigenstepper_targetSurface_and_onSurface,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    using ActionListType      = ActionList<>;
    using AbortConditionsType = AbortList<>;
    using PropagateState      = MultiEigenPropagatorType::
        State<PropagatorOptions<ActionListType, AbortConditionsType>>;

    // setup propagation options
    PropagatorOptions<ActionListType, AbortConditionsType> options(tgContext,
                                                                   mfContext);

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    // the 1st component
    Vector3D mom_1(px, py, pz);
    // the 2nd component
    Vector3D mom_2(px + 0.5 * px, py, pz);
    // flip the mom(dir) to make a normal component that can not reach the
    // target surface
    Vector3D                         mom_3 = -1 * mom_1;
    CurvilinearParameters            start1(nullptr, pos, mom_1, q);
    CurvilinearParameters            start2(nullptr, pos, mom_2, q);
    CurvilinearParameters            start3(nullptr, pos, mom_3, q);
    PropagateState                   prop_state(start1, options);
    VoidIntersectionCorrector        corrector;
    Acts::NavigationOptions<Surface> opt(forward, true);

    // append a component to the end of the column
    // the component has different mom(dir)
    prop_state.stepping.stateCol.push_back(std::make_tuple(
        MultiEigenStepperType::SingleStateType(
            tgContext, mfContext, start2, forward, options.maxStepSize),
        1.,
        StateStatus::Free));
    // test the normalize function, each of components is weighted 0.5
    mcstepper.normalizeComponents(prop_state.stepping);
    for (auto& tuple_state : prop_state.stepping.stateCol) {
      auto weight = std::get<1>(tuple_state);
      CHECK_CLOSE_REL(weight, 0.5, 0.0001);
    }

    /// Four steps to pass through surface a,b,c
    /// 1. propagate to surface-a, test if all components calculate their
    /// stepSize alone,
    ///    so that they reach the surface in a right way;
    ///    the surfaceReached method test if they are on surface, if so then set
    ///    them free.
    /// 2. propagate to surface-b, mandatorily change the stepSize of 1st
    /// component to make it fail to reach the surface.
    /// 3. propagate to surface-b again, virtually only propagate the 1st
    /// component since 2nd has been Locked on the surface-b.
    /// 4. propagate to surface-c with adding a wired component with small
    /// weight and a flip direction,
    ///    this component is tested to be killed since it will never reach the
    ///    target surface.

    /// First, propagate to surface-a
    // update the stepSize of each components to surface-a
    mcstepper.targetSurface(
        prop_state.stepping, aSurface.get(), opt, corrector);
    // step to surface-a
    mcstepper.step(prop_state);

    bool isAllOnSurface
        = mcstepper.surfaceReached(prop_state.stepping, aSurface.get());
    // all components should be on surface-a, and then release as free
    // components
    BOOST_CHECK(isAllOnSurface);

    // the two components should both be arrived at the a-surface
    // the status now should be {Free, Free}( {0, 0} )
    // test if stepper.surfaceReached() can get the right result
    // test if mcstepper.surfaceReached() lock the right component at the
    // current surface
    std::array<int, 2> reference_status_check  = {0, 0};
    std::array<int, 2> reference_surface_check = {1, 1};
    std::array<int, 2> status_check;
    std::array<int, 2> surface_check;
    unsigned int id = 0;
    for (auto& tuple_state : prop_state.stepping.stateCol) {
      bool isOnSurface
          = stepper.surfaceReached(std::get<0>(tuple_state), aSurface.get());
      surface_check[id] = isOnSurface;
      status_check[id]  = int(std::get<2>(tuple_state));
      id++;
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_status_check.begin(),
                                  reference_status_check.end(),
                                  status_check.begin(),
                                  status_check.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_surface_check.begin(),
                                  reference_surface_check.end(),
                                  surface_check.begin(),
                                  surface_check.end());

    /// Secondly, propagate to surface-b, this time force the 1st component
    /// stepSize to avoid it reaching that surface
    // update the stepSize of each components to surface-b
    mcstepper.targetSurface(
        prop_state.stepping, bSurface.get(), opt, corrector);
    // manually minimize the stepSize of 1st component so that it can not reach
    // the surface
    auto& firstComponent = std::get<0>(*prop_state.stepping.stateCol.begin());
    firstComponent.stepSize.update(
        firstComponent.stepSize.value(cstep::actor) / 10, cstep::actor);
    // step to surface-b
    mcstepper.step(prop_state);

    isAllOnSurface
        = mcstepper.surfaceReached(prop_state.stepping, bSurface.get());
    // the 1st component does not reach, so it should return false
    BOOST_CHECK(!isAllOnSurface);

    // the 1st component should not reach the surface, and the 2nd component
    // should reach the surface at this time
    // the status now should be {Free, Locked}
    // test if stepper.surfaceReached() can get the right result
    // test if mcstepper.surfaceReached() lock the right component at the
    // current surface
    reference_status_check  = {0, 1};
    reference_surface_check = {0, 1};
    id                      = 0;
    for (auto& tuple_state : prop_state.stepping.stateCol) {
      bool isOnSurface
          = stepper.surfaceReached(std::get<0>(tuple_state), bSurface.get());
      surface_check[id] = isOnSurface;
      status_check[id]  = int(std::get<2>(tuple_state));
      id++;
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_status_check.begin(),
                                  reference_status_check.end(),
                                  status_check.begin(),
                                  status_check.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_surface_check.begin(),
                                  reference_surface_check.end(),
                                  surface_check.begin(),
                                  surface_check.end());

    /// Then, propagate to surface-b again
    // update the stepSize of each components to surface-b
    mcstepper.targetSurface(
        prop_state.stepping, bSurface.get(), opt, corrector);
    // step to surface-b again
    mcstepper.step(prop_state);
    isAllOnSurface
        = mcstepper.surfaceReached(prop_state.stepping, bSurface.get());
    BOOST_CHECK(isAllOnSurface);
    reference_status_check  = {0, 0};
    reference_surface_check = {1, 1};
    id                      = 0;
    for (auto& tuple_state : prop_state.stepping.stateCol) {
      bool isOnSurface
          = stepper.surfaceReached(std::get<0>(tuple_state), bSurface.get());
      surface_check[id] = isOnSurface;
      status_check[id]  = int(std::get<2>(tuple_state));
      id++;
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_status_check.begin(),
                                  reference_status_check.end(),
                                  status_check.begin(),
                                  status_check.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_surface_check.begin(),
                                  reference_surface_check.end(),
                                  surface_check.begin(),
                                  surface_check.end());

    // add a component of which direction is opposite but with very small weight
    // - so that not influenced the direction stream
    prop_state.stepping.stateCol.push_back(std::make_tuple(
        MultiEigenStepperType::SingleStateType(
            tgContext, mfContext, start3, forward, options.maxStepSize),
        10e-6,
        StateStatus::Free));
    mcstepper.normalizeComponents(prop_state.stepping);

    /// Finally, propagate to surface-c, this time the new added component can
    /// never reach surface-c, so it will be killed
    // update the stepSize of each components to surface-c, meanwhile the wierd
    // component is killed
    mcstepper.targetSurface(
        prop_state.stepping, cSurface.get(), opt, corrector);
    BOOST_CHECK_EQUAL(prop_state.stepping.stateCol.size(), 2);
    // step to surface-c
    mcstepper.step(prop_state);

    isAllOnSurface
        = mcstepper.surfaceReached(prop_state.stepping, cSurface.get());
    BOOST_CHECK(isAllOnSurface);

    reference_status_check  = {0, 0};
    reference_surface_check = {1, 1};
    id                      = 0;
    for (auto& tuple_state : prop_state.stepping.stateCol) {
      // test if the remaining components both own their previous weights 0.5
      CHECK_CLOSE_REL(std::get<1>(tuple_state), 0.5, 0.0001);
      bool isOnSurface
          = stepper.surfaceReached(std::get<0>(tuple_state), cSurface.get());
      surface_check[id] = isOnSurface;
      status_check[id]  = int(std::get<2>(tuple_state));
      id++;
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_status_check.begin(),
                                  reference_status_check.end(),
                                  status_check.begin(),
                                  status_check.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(reference_surface_check.begin(),
                                  reference_surface_check.end(),
                                  surface_check.begin(),
                                  surface_check.end());
  }

}  // namespace Test
}  // namespace Acts
