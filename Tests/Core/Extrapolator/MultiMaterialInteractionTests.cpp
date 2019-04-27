// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE MultiMaterialInteraction Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiEigenStepper.hpp"
#include "Acts/Extrapolator/MultiMaterialInteractor.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Utilities/MagneticFieldContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

namespace Acts {
namespace Test {

  using DebugOutput = detail::DebugOutputActor;

  std::normal_distribution<double> gauss(0., 1.);
  std::default_random_engine       generator(42);

  bool debugMode = false;

  // Create a test context
  GeometryContext      tgContext  = GeometryContext();
  MagneticFieldContext mfContext  = MagneticFieldContext();
  CalibrationContext   calContext = CalibrationContext();

  ///
  /// @brief Unit test for multi material interactor
  /// This tests two things:
  /// 1.The interacion points are the same at each surfaces when propagate,
  /// since multimaterial interaction splits into 2 same components currently.
  /// 2.With the cubic box of 6 surfaces, the final number of components is
  /// equal to 2^6.
  BOOST_AUTO_TEST_CASE(multi_material_interactor)
  {
    // Build detector
    CubicTrackingGeometry cGeometry(tgContext);
    auto                  detector = cGeometry();
    // Build navigator
    Navigator mNavigator(detector);
    mNavigator.resolvePassive   = false;
    mNavigator.resolveMaterial  = true;
    mNavigator.resolveSensitive = true;

    ConstantBField bField(Vector3D(0., 0., 0.));
    using MultiEigenStepperType = MultiEigenStepper<ConstantBField>;
    MultiEigenStepperType multi_stepper(bField);
    using MultiPropagator = Propagator<MultiEigenStepperType, Navigator>;
    MultiPropagator multiPropagator(multi_stepper, mNavigator);
    //
    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 1000. * units::_um, 0., 0., 0., 0., 0., 1000. * units::_um, 0., 0.,
        0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.01;

    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    Vector3D rPos(-3. * units::_m,
                  10. * units::_um * gauss(generator),
                  100. * units::_um * gauss(generator));
    Vector3D rMom(1. * units::_GeV,
                  0.025 * units::_GeV * gauss(generator),
                  0.025 * units::_GeV * gauss(generator));

    SingleCurvilinearTrackParameters<ChargedPolicy> rStart(
        std::move(covPtr), rPos, rMom, 1.);

    PropagatorOptions<ActionList<DebugOutput, MultiMaterialInteractor>,
                      AbortList<detail::EndOfWorldReached>>
        rOptions(tgContext, mfContext);
    rOptions.debug = debugMode;

    auto result = multiPropagator.propagate(rStart, rOptions).value();
    auto numOfComponents
        = result.template get<MultiMaterialInteractor::result_type>()
              .numComponents;
    if (debugMode) {
      const auto debugString
          = result.template get<DebugOutput::result_type>().debugString;
      std::cout << ">>>> Measurement creation: " << std::endl;
      std::cout << debugString;
      std::cout << " In the Propagator, 1 component finally splits into "
                << numOfComponents << " components.";
    }

    // Since The emptyInteractor 1 component split into 2 same components
    // Test if the number of components split into 64 in the interactions of 6
    // surfaces of Cubic Box
    BOOST_CHECK(numOfComponents == 64);

    // Test at each surface if all material interactions recorded are the same
    // because the component split don't change anything, just copy component at
    // current status
    const auto& material_interactions_result
        = result.template get<MultiMaterialInteractor::result_type>()
              .multiMaterialInteractions;
    for (const auto& materialInteractionPair : material_interactions_result) {
      const auto& materialInteractionVec = materialInteractionPair.second;
      BOOST_CHECK(
          std::all_of(materialInteractionVec.begin() + 1,
                      materialInteractionVec.end(),
                      [&](const Acts::InteractionPointVec::value_type& r) {
                        return r == materialInteractionVec.front();
                      }));
    }
  }
}
}
