// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief The Material interaction struct
/// It records the surface  and the passed material
/// This is only nessecary recorded when configured
struct MaterialInteraction {
  /// The material surface
  const Surface* surface = nullptr;

  double sigmaPhi2 = 0.;    ///< applied material effect: sigma(phi)^2
  double sigmaTheta2 = 0.;  ///< applied material effect: sigma(theta)^2
  double deltaP = 0.;       ///< applied material effect: dela(p)
  double sigmaQoP2 = 0.;    ///< applied material effect: sigma(qop)^2

  /// The position information of the material hit
  Vector3D position = Vector3D(0., 0., 0);
  /// The direction information of the material hit
  Vector3D direction = Vector3D(0., 0., 0);
  /// The calculated path & applied path correction factor
  double pathCorrection = 1.;
  /// The (passsed) material properties
  /// it is the material and the actual (corrected) path length
  MaterialProperties materialProperties = MaterialProperties();
};

/// The Material interactor struct
///
/// This is a plugin to the Propagator that
/// performs material interaction on the currentSurface
/// of the Propagagor state
struct MaterialInteractor {
  // Configuration for this MaterialInteractor
  detail::PointwiseMaterialInteraction interaction;

  /// multiple scattering switch on/off
  bool multipleScattering = true;

  /// Energy loss switch on/off
  bool energyLoss = true;

  /// Record material in detail
  bool recordInteractions = false;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state which is
  /// created for every propagation/extrapolation step
  struct this_result {
    // The accumulated materialInX0
    double materialInX0 = 0.;
    /// The accumulated materialInL0
    double materialInL0 = 0.;
    /// This one is only filled when recordInteractions is switched on
    std::vector<MaterialInteraction> materialInteractions;
  };

  using result_type = this_result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    using namespace Acts::UnitLiterals;

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    // If switched off, then return - alows run-time configuration
    if (!multipleScattering && !energyLoss && !recordInteractions) {
      return;
    }

    // A current surface has been already assigned by the navigator
    // check for material
    if (state.navigation.currentSurface &&
        state.navigation.currentSurface->surfaceMaterial()) {
      MaterialProperties mProperties =
          interaction.evaluateMaterialProperties(state, stepper);
      // Material properties (non-zero) have been found for this configuration
      if (mProperties) {
        // more debugging output to the screen
        debugLog(state, [&] {
          return std::string("Material properties found for this surface.");
        });

        // To integrate process noise, we need to transport
        // the covariance to the current position in space
        // the 'true' indicates re-initializaiton of the further transport
        if (state.stepping.covTransport) {
          stepper.covarianceTransport(state.stepping, true);
        }

        // Calculate the path correction
        double pCorrection = state.navigation.currentSurface->pathCorrection(
            state.geoContext, stepper.position(state.stepping),
            stepper.direction(state.stepping));

        // Scale the material properties
        mProperties *= pCorrection;

        // Perform the material interactions and collect Var[phi], Var[theta],
        // DeltaP, Var[QoP]
        std::array<double, 4> sigma = interaction.evaluateMaterialInteraction(
            state, stepper, mProperties, multipleScattering, energyLoss);

        if (multipleScattering) {
          applyMS(state, sigma[0], sigma[1]);
        }
        if (energyLoss) {
          applyEnergyLoss(state, stepper, sigma[2], sigma[3]);
        }

        // This doesn't cost anything - do it regardless
        result.materialInX0 += mProperties.thicknessInX0();
        result.materialInL0 += mProperties.thicknessInL0();

        // Record the material interaction if configured to do so
        if (recordInteractions) {
          // Create the material interaction class, in case we record afterwards
          MaterialInteraction mInteraction;
          mInteraction.surface = state.navigation.currentSurface;
          // Record the material interaction
          mInteraction.sigmaPhi2 = sigma[0];
          mInteraction.sigmaTheta2 = sigma[1];
          mInteraction.deltaP = sigma[2];
          mInteraction.sigmaQoP2 = sigma[3];

          mInteraction.position = stepper.position(state.stepping);
          mInteraction.direction = stepper.direction(state.stepping);
          mInteraction.materialProperties = mProperties;
          mInteraction.pathCorrection = pCorrection;
          result.materialInteractions.push_back(std::move(mInteraction));
        }
      }
    }
  }

  /// Pure observer interface
  /// This does not apply to the surface collector
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /*state*/) const {}

 private:
  /// @brief This function applies the changes to the state due to multiple
  /// scattering
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  ///
  /// @param [in, out] state is the mutable propagator state object
  /// @param [in] sigmaPhi2 The variance in phi by multiple scattering
  /// @param [in] sigmaTheta2 The variance in theta by multiple scattering
  template <typename propagator_state_t>
  void applyMS(propagator_state_t& state, double sigmaPhi2,
               double sigmaTheta2) const {
    if (state.stepping.navDir == forward) {
      // Just add the multiple scattering component
      state.stepping.cov(ePHI, ePHI) += state.stepping.navDir * sigmaPhi2;
      state.stepping.cov(eTHETA, eTHETA) += state.stepping.navDir * sigmaTheta2;
    } else {
      // We check if the covariance stays positive
      const double sEphi = state.stepping.cov(ePHI, ePHI);
      const double sEtheta = state.stepping.cov(eTHETA, eTHETA);
      if (sEphi > sigmaPhi2 && sEtheta > sigmaTheta2) {
        // Noise removal is not applied if covariance would fall below 0
        state.stepping.cov(ePHI, ePHI) -= sigmaPhi2;
        state.stepping.cov(eTHETA, eTHETA) -= sigmaTheta2;
      }
    }
  }

  /// @brief This function applies the changes to the state due to energy loss
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param [in, out] state is the mutable propagator state object
  /// @param [in] stepper The stepper in use
  /// @param [in] deltaP The momentum change
  /// @param [in] sigmaQoP2 The variance in QoP by energy loss
  template <typename propagator_state_t, typename stepper_t>
  void applyEnergyLoss(propagator_state_t& state, const stepper_t& stepper,
                       double deltaP, double sigmaQoP2) const {
    // Update the state/momentum
    stepper.update(state.stepping, stepper.position(state.stepping),
                   stepper.direction(state.stepping),
                   std::copysign(stepper.momentum(state.stepping) - deltaP,
                                 stepper.momentum(state.stepping)),
                   stepper.time(state.stepping));

    // Good in any case for positive direction
    if (state.stepping.navDir == forward) {
      state.stepping.cov(eQOP, eQOP) += state.stepping.navDir * sigmaQoP2;
    } else {
      // Check that covariance entry doesn't become negative
      const double sEqop = state.stepping.cov(eQOP, eQOP);
      if (sEqop > sigmaQoP2) {
        state.stepping.cov(eQOP, eQOP) += state.stepping.navDir * sigmaQoP2;
      }
    }
  }

  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param state the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a streamable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(state.options.debugPfxWidth);
      dstream << "material interaction"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};

/// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::result_type;

/// And recorded material track
/// - this is start:  position, start momentum
///   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;

}  // end of namespace Acts