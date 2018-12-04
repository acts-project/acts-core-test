// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Stepper Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <fstream>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

  using cstep = detail::ConstrainedStep;

  ///
  /// @brief Aborter for the case that a particle leaves the detector or reaches
  /// a custom made threshold.
  ///
  struct EndOfWorld
  {
    /// Maximum value in x-direction of the detector
    double maxX = 1. * units::_m;

    /// @brief Constructor
    EndOfWorld() = default;

    /// @brief Main call operator for the abort operation
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    /// @return Boolean statement if the particle is still in the detector
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (std::abs(state.stepping.position().x()) >= maxX
          || std::abs(state.stepping.position().y()) >= 0.5 * units::_m
          || std::abs(state.stepping.position().z()) >= 0.5 * units::_m)
        return true;
      return false;
    }
  };

  ///
  /// @brief Data collector while propagation
  ///
  struct StepCollector
  {

    ///
    /// @brief Data container for result analysis
    ///
    struct this_result
    {
      // Position of the propagator after each step
      std::vector<Vector3D> position;
      // Momentum of the propagator after each step
      std::vector<Vector3D> momentum;
    };

    using result_type = this_result;

    /// @brief Main call operator for the action list. It stores the data for
    /// analysis afterwards
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @param [in] state State of the propagator
    /// @param [out] result Struct which is filled with the data
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      result.position.push_back(state.stepping.position());
      result.momentum.push_back(state.stepping.momentum());
    }
  };

  /// @brief This function tests the EigenStepper with the DefaultExtension and
  /// the DenseEnvironmentExtension. The focus of this tests lies in the
  /// choosing of the right extension for the individual use case. This is
  /// performed with three different detectors:
  /// a) Pure vaccuum -> DefaultExtension needs to act
  /// b) Pure Be -> DenseEnvironmentExtension needs to act
  /// c) Vacuum - Be - Vacuum -> Both should act and switch during the
  /// propagation

  // Test case a). The DenseEnvironmentExtension should state that it is not
  // valid in this case.
  BOOST_AUTO_TEST_CASE(step_extension_vacuum_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0.5 * units::_m, 0., 0.};
    vConf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {0.5 * units::_m, 0., 0.};
    conf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> vacuum = tgb.trackingGeometry();

    // Build navigator
    Navigator naviVac(vacuum);
    naviVac.resolvePassive   = true;
    naviVac.resolveMaterial  = true;
    naviVac.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList     = aList;
    propOpts.stopConditions = abortList;
    propOpts.maxSteps       = 1e6;
    propOpts.maxStepSize    = 0.5 * units::_m;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviVac);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Check that the propagation happend without interactions
    for (const auto& pos : stepResult.position) {
      BOOST_TEST(pos.y() == 0.);
      BOOST_TEST(pos.z() == 0.);
      if (pos == stepResult.position.back())
        BOOST_TEST(pos.x() == 1. * units::_m);
    }
    for (const auto& mom : stepResult.momentum) {
      BOOST_TEST(mom.x() == 1. * units::_GeV);
      BOOST_TEST(mom.y() == 0.);
      BOOST_TEST(mom.z() == 0.);
    }

    // Rebuild and check the choice of extension
    ActionList<StepCollector> aListDef;

    // Set options for propagator
    PropagatorOptions<ActionList<StepCollector>, AbortList<EndOfWorld>>
        propOptsDef;
    propOptsDef.actionList     = aListDef;
    propOptsDef.stopConditions = abortList;
    propOptsDef.maxSteps       = 1e6;
    propOptsDef.maxStepSize    = 0.5 * units::_m;

    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension>>
        esDef(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension>>,
               Navigator>
        propDef(esDef, naviVac);

    // Launch and collect results
    const auto& resultDef = propDef.propagate(sbtp, propOptsDef);
    const StepCollector::this_result& stepResultDef
        = resultDef.get<typename StepCollector::result_type>();

    // Check that the right extension was chosen
    // If chosen correctly, the number of elements should be identical
    BOOST_TEST(stepResult.position.size() == stepResultDef.position.size());
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      BOOST_TEST(stepResult.position[i] == stepResultDef.position[i]);
    }
    BOOST_TEST(stepResult.momentum.size() == stepResultDef.momentum.size());
    for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
      BOOST_TEST(stepResult.momentum[i] == stepResultDef.momentum[i]);
    }
  }
  // Test case b). The DefaultExtension should state that it is invalid here.
  BOOST_AUTO_TEST_CASE(step_extension_material_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0.5 * units::_m, 0., 0.};
    vConf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf.material = std::make_shared<const Material>(
        Material(352.8, 394.133, 9.012, 4., 1.848e-3));
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {0.5 * units::_m, 0., 0.};
    conf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> material = tgb.trackingGeometry();

    // Build navigator
    Navigator naviMat(material);
    naviMat.resolvePassive   = true;
    naviMat.resolveMaterial  = true;
    naviMat.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(5. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList     = aList;
    propOpts.stopConditions = abortList;
    propOpts.maxSteps       = 1e6;
    propOpts.maxStepSize    = 0.5 * units::_m;
    propOpts.debug          = true;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviMat);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Check that there occured interaction
    for (const auto& pos : stepResult.position) {
      BOOST_TEST(pos.y() == 0.);
      BOOST_TEST(pos.z() == 0.);
      if (pos == stepResult.position.front()) {
        BOOST_TEST(pos.x() == 0.);
      } else {
        BOOST_TEST(pos.x() != 0.);
      }
    }
    for (const auto& mom : stepResult.momentum) {
      BOOST_TEST(mom.y() == 0.);
      BOOST_TEST(mom.z() == 0.);
      if (mom == stepResult.momentum.front()) {
        BOOST_TEST(mom.x() == 5. * units::_GeV);
      } else {
        BOOST_TEST(mom.x() < 5. * units::_GeV);
      }
    }

    // Rebuild and check the choice of extension
    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOptsDense;
    propOptsDense.actionList     = aList;
    propOptsDense.stopConditions = abortList;
    propOptsDense.maxSteps       = 1e6;
    propOptsDense.maxStepSize    = 0.5 * units::_m;
    propOptsDense.debug          = true;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DenseEnvironmentExtension>>
        esDense(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DenseEnvironmentExtension>>,
               Navigator>
        propDense(esDense, naviMat);

    // Launch and collect results
    const auto& resultDense = propDense.propagate(sbtp, propOptsDense);
    const StepCollector::this_result& stepResultDense
        = resultDense.get<typename StepCollector::result_type>();

    // Check that the right extension was chosen
    // If chosen correctly, the number of elements should be identical
    BOOST_TEST(stepResult.position.size() == stepResultDense.position.size());
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      BOOST_TEST(stepResult.position[i] == stepResultDense.position[i]);
    }
    BOOST_TEST(stepResult.momentum.size() == stepResultDense.momentum.size());
    for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
      BOOST_TEST(stepResult.momentum[i] == stepResultDense.momentum[i]);
    }

    ////////////////////////////////////////////////////////////////////

    // Re-launch the configuration with magnetic field
    bField.setField(0., 1. * units::_T, 0.);
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        esB(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        propB(esB, naviMat);

    const auto& resultB = propB.propagate(sbtp, propOptsDense);
    const StepCollector::this_result& stepResultB
        = resultB.get<typename StepCollector::result_type>();

    // Check that there occured interaction
    for (const auto& pos : stepResultB.position) {
      if (pos == stepResultB.position.front()) {
        BOOST_TEST(pos.x() == 0.);
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() == 0.);
      } else {
        BOOST_TEST(pos.x() != 0.);
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() != 0.);
      }
    }
    for (const auto& mom : stepResultB.momentum) {
      if (mom == stepResultB.momentum.front()) {
        BOOST_TEST(mom.x() == 5. * units::_GeV);
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() == 0.);
      } else {
        BOOST_TEST(mom.x() != 5. * units::_GeV);
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() != 0.);
      }
    }
  }
  // Test case c). Both should be involved in their part of the detector
  BOOST_AUTO_TEST_CASE(step_extension_vacmatvac_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConfVac1;
    vConfVac1.position = {0.5 * units::_m, 0., 0.};
    vConfVac1.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    CuboidVolumeBuilder::VolumeConfig vConfMat;
    vConfMat.position = {1.5 * units::_m, 0., 0.};
    vConfMat.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConfMat.material = std::make_shared<const Material>(
        Material(352.8, 394.133, 9.012, 4., 1.848e-3));
    CuboidVolumeBuilder::VolumeConfig vConfVac2;
    vConfVac2.position = {2.5 * units::_m, 0., 0.};
    vConfVac2.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg = {vConfVac1, vConfMat, vConfVac2};
    conf.position  = {1.5 * units::_m, 0., 0.};
    conf.length    = {3. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> det = tgb.trackingGeometry();

    // Build navigator
    Navigator naviDet(det);
    naviDet.resolvePassive   = true;
    naviDet.resolveMaterial  = true;
    naviDet.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(5. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;
    abortList.get<EndOfWorld>().maxX = 3. * units::_m;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList     = aList;
    propOpts.stopConditions = abortList;
    propOpts.maxSteps       = 1e6;
    propOpts.maxStepSize    = 0.5 * units::_m;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 1. * units::_T, 0.));
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviDet);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Manually set the extensions for each step and propagate through each
    // volume by propagation to the boundaries
    // Collect boundaries
    std::vector<Surface const*> surs;
    std::vector<std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>>
        boundaries = det->lowestTrackingVolume({0.5 * units::_m, 0., 0.})
                         ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 1. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }
    boundaries = det->lowestTrackingVolume({1.5 * units::_m, 0., 0.})
                     ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 2. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }
    boundaries = det->lowestTrackingVolume({2.5 * units::_m, 0., 0.})
                     ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 3. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }

    // Build launcher through vacuum
    // Set options for propagator
    ActionList<StepCollector> aListDef;

    PropagatorOptions<ActionList<StepCollector>, AbortList<EndOfWorld>>
        propOptsDef;
    abortList.get<EndOfWorld>().maxX = 1. * units::_m;
    propOptsDef.actionList           = aListDef;
    propOptsDef.stopConditions       = abortList;
    propOptsDef.maxSteps             = 1e6;
    propOptsDef.maxStepSize          = 0.5 * units::_m;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension>>
        esDef(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension>>,
               Navigator>
        propDef(esDef, naviDet);

    // Launch and collect results
    const auto& resultDef = propDef.propagate(sbtp, *(surs[0]), propOptsDef);
    const StepCollector::this_result& stepResultDef
        = resultDef.get<typename StepCollector::result_type>();

    // Check the exit situation of the first volume
    std::pair<Vector3D, Vector3D> endParams, endParamsControl;
    for (unsigned int i = 0; i < stepResultDef.position.size(); i++) {
      if (1. * units::_m - stepResultDef.position[i].x() < 1e-4) {
        endParams = std::make_pair(stepResultDef.position[i],
                                   stepResultDef.momentum[i]);
        break;
      }
    }
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (1. * units::_m - stepResult.position[i].x() < 1e-4) {
        endParamsControl
            = std::make_pair(stepResult.position[i], stepResult.momentum[i]);
        break;
      }
    }

    BOOST_TEST(endParams.first.x() == endParamsControl.first.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.y() == endParamsControl.first.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.z() == endParamsControl.first.z(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.x() == endParamsControl.second.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.y() == endParamsControl.second.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.z() == endParamsControl.second.z(),
               tt::tolerance(1e-5));

    BOOST_TEST(endParams.first.x() == endParamsControl.first.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.y() == endParamsControl.first.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.z() == endParamsControl.first.z(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.x() == endParamsControl.second.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.y() == endParamsControl.second.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.z() == endParamsControl.second.z(),
               tt::tolerance(1e-5));

    // Build launcher through material
    // Set initial parameters for the particle track by using the result of the
    // first volume
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(
        ActsSymMatrixD<5>::Identity());
    startParams = endParams.first;
    startMom    = endParams.second;
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtpPiecewise(
        std::move(covPtr), startParams, startMom, 1.);

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOptsDense;
    abortList.get<EndOfWorld>().maxX = 2. * units::_m;
    propOptsDense.actionList         = aList;
    propOptsDense.stopConditions     = abortList;
    propOptsDense.maxSteps           = 1e6;
    propOptsDense.maxStepSize        = 0.5 * units::_m;
    propOptsDense.tolerance          = 1e-8;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DenseEnvironmentExtension>>
        esDense(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DenseEnvironmentExtension>>,
               Navigator>
        propDense(esDense, naviDet);

    // Launch and collect results
    const auto& resultDense
        = propDense.propagate(sbtpPiecewise, *(surs[1]), propOptsDense);
    const StepCollector::this_result& stepResultDense
        = resultDense.get<typename StepCollector::result_type>();

    // Check the exit situation of the second volume
    for (unsigned int i = 0; i < stepResultDense.position.size(); i++) {
      if (2. * units::_m - stepResultDense.position[i].x() < 1e-4) {
        endParams = std::make_pair(stepResultDense.position[i],
                                   stepResultDense.momentum[i]);
        break;
      }
    }
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (2. * units::_m - stepResult.position[i].x() < 1e-4) {
        endParamsControl
            = std::make_pair(stepResult.position[i], stepResult.momentum[i]);
        break;
      }
    }

    BOOST_TEST(endParams.first.x() == endParamsControl.first.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.y() == endParamsControl.first.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.first.z() == endParamsControl.first.z(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.x() == endParamsControl.second.x(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.y() == endParamsControl.second.y(),
               tt::tolerance(1e-5));
    BOOST_TEST(endParams.second.z() == endParamsControl.second.z(),
               tt::tolerance(1e-5));
  }

  // Test case a). The DenseEnvironmentExtension should state that it is not
  // valid in this case.
  BOOST_AUTO_TEST_CASE(step_extension_trackercalomdt_test)
  {
    double             rotationAngle = M_PI * 0.5;
    Vector3D           xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D           yPos(0., 1., 0.);
    Vector3D           zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    MaterialProperties matProp(
        352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);

    CuboidVolumeBuilder                cvb;
    CuboidVolumeBuilder::SurfaceConfig sConf1;
    sConf1.position        = Vector3D(0.3 * units::_m, 0., 0.);
    sConf1.rotation.col(0) = xPos;
    sConf1.rotation.col(1) = yPos;
    sConf1.rotation.col(2) = zPos;
    sConf1.rBounds         = std::make_shared<const RectangleBounds>(
        RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    sConf1.surMat = std::shared_ptr<const SurfaceMaterial>(
        new HomogeneousSurfaceMaterial(matProp));
    sConf1.thickness = 1. * units::_mm;
    CuboidVolumeBuilder::LayerConfig lConf1;
    lConf1.surfaceCfg = sConf1;

    CuboidVolumeBuilder::SurfaceConfig sConf2;
    sConf2.position        = Vector3D(0.6 * units::_m, 0., 0.);
    sConf2.rotation.col(0) = xPos;
    sConf2.rotation.col(1) = yPos;
    sConf2.rotation.col(2) = zPos;
    sConf2.rBounds         = std::make_shared<const RectangleBounds>(
        RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    sConf2.surMat = std::shared_ptr<const SurfaceMaterial>(
        new HomogeneousSurfaceMaterial(matProp));
    sConf2.thickness = 1. * units::_mm;
    CuboidVolumeBuilder::LayerConfig lConf2;
    lConf2.surfaceCfg = sConf2;

    CuboidVolumeBuilder::VolumeConfig muConf1;
    muConf1.position = {2.3 * units::_m, 0., 0.};
    muConf1.length   = {20. * units::_cm, 20. * units::_cm, 20. * units::_cm};
    muConf1.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    muConf1.name = "MDT1";
    CuboidVolumeBuilder::VolumeConfig muConf2;
    muConf2.position = {2.7 * units::_m, 0., 0.};
    muConf2.length   = {20. * units::_cm, 20. * units::_cm, 20. * units::_cm};
    muConf2.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    muConf2.name = "MDT2";

    CuboidVolumeBuilder::VolumeConfig vConf1;
    vConf1.position = {0.5 * units::_m, 0., 0.};
    vConf1.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf1.layerCfg = {lConf1, lConf2};
    vConf1.name     = "Tracker";
    CuboidVolumeBuilder::VolumeConfig vConf2;
    vConf2.position = {1.5 * units::_m, 0., 0.};
    vConf2.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf2.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    vConf2.name = "Calorimeter";
    CuboidVolumeBuilder::VolumeConfig vConf3;
    vConf3.position  = {2.5 * units::_m, 0., 0.};
    vConf3.length    = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf3.volumeCfg = {muConf1, muConf2};
    vConf3.name      = "Muon system";
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg = {vConf1, vConf2, vConf3};
    conf.position  = {1.5 * units::_m, 0., 0.};
    conf.length    = {3. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> detector = tgb.trackingGeometry();

    // Build navigator
    Navigator naviVac(detector);
    naviVac.resolvePassive   = true;
    naviVac.resolveMaterial  = true;
    naviVac.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector, MaterialInteractor>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.stopConditions.get<EndOfWorld>().maxX = 3. * units::_m;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField,
                 VoidCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviVac);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Test that momentum changes only occured at the right detector parts
    double lastMomentum = stepResult.momentum[0].x();
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      // Test for changes
      if ((stepResult.position[i].x() > 0.3 * units::_m
           && stepResult.position[i].x() < 0.6 * units::_m)
          || (stepResult.position[i].x() > 0.6 * units::_m
              && stepResult.position[i].x() <= 1. * units::_m)
          || (stepResult.position[i].x() > 1. * units::_m
              && stepResult.position[i].x() <= 2. * units::_m)
          || (stepResult.position[i].x() > 2.2 * units::_m
              && stepResult.position[i].x() <= 2.4 * units::_m)
          || (stepResult.position[i].x() > 2.6 * units::_m
              && stepResult.position[i].x() <= 2.8 * units::_m)) {
        BOOST_TEST(stepResult.momentum[i].x() <= lastMomentum);
        lastMomentum = stepResult.momentum[i].x();
      } else
      // Test the absence of momentum loss
      {
        if (stepResult.position[i].x() < 0.3 * units::_m
            || (stepResult.position[i].x() > 2. * units::_m
                && stepResult.position[i].x() <= 2.2 * units::_m)
            || (stepResult.position[i].x() > 2.4 * units::_m
                && stepResult.position[i].x() <= 2.6 * units::_m)
            || (stepResult.position[i].x() > 2.8 * units::_m
                && stepResult.position[i].x() <= 3. * units::_m)) {
          BOOST_TEST(stepResult.momentum[i].x() == lastMomentum);
        }
      }
    }
  }
}  // namespace Test
}  // namespace Acts
