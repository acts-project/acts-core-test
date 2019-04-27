// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Extrapolator/detail/MaterialEffect.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/detail/InteractionFormulas.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @brief The Material interaction struct
/// It records the surface and the passed material
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

  bool
  operator==(const MaterialInteraction& others) const
  {
    if (fabs((this->direction - others.direction).norm()) > 1e-10
        || fabs((this->position - others.position).norm()) > 1e-10
        || this->surface != others.surface) {
      return false;
    }
    return true;
  }
};

/// The Material interactor struct
///
/// This is a plugin to the Propagator that
/// performs material interaction on the currentSurface
/// of the Propagagor state
struct MaterialInteractor {
  // Configuration for this MaterialInteractor

  /// multiple scattering switch on/off
  bool multipleScattering = true;
  /// The scattering formula struct
  detail::HighlandScattering process_scattering;
  /// apply delta theta/phi on the covariance matrix
  detail::MultipleScattering multiscattering;

  /// Energy loss switch on/off
  bool energyLoss = true;
  /// The energy loss formula struct
  detail::IonisationLoss process_ionisationloss;
  /// apply dE and delta covariance on the matrix
  detail::EnergyLoss energyloss;

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
      // Let's set the pre/full/post update stage
      MaterialUpdateStage mStage = fullUpdate;
      // We are at the start surface
      if (state.navigation.startSurface == state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on start surface: post-update mode.");
        });
        mStage = postUpdate;
        // Or is it the target surface ?
      } else if (state.navigation.targetSurface ==
                 state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on target surface: pre-update mode");
        });
        mStage = preUpdate;
      } else {
        debugLog(state, [&] {
          return std::string("Update while pass through: full mode.");
        });
      }

      // Get the surface material & properties from them and continue if you
      // found some
      const ISurfaceMaterial* sMaterial =
          state.navigation.currentSurface->surfaceMaterial();
      MaterialProperties mProperties = sMaterial->materialProperties(
          stepper.position(state.stepping), state.stepping.navDir, mStage);
      // Material properties (non-zero) have been found for this configuration
      if (mProperties) {
        // more debugging output to the screen
        debugLog(state, [&] {
          return std::string("Material properties found for this surface.");
        });

        // Create the material interaction class, in case we record afterwards
        MaterialInteraction mInteraction;
        mInteraction.surface = state.navigation.currentSurface;

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

        // The momentum at current position
        const double p = stepper.momentum(state.stepping);
        const double m = state.options.mass;
        const double E = std::sqrt(p * p + m * m);
        const double lbeta = p / E;

        // Apply the multiple scattering
        // - only when you do covariance transport
        if (multipleScattering && state.stepping.covTransport) {
          // Thickness in X0 from without path correction
          double tInX0 = mProperties.thicknessInX0();
          // Retrieve the scattering contribution
          double sigmaScat = process_scattering(p, lbeta, tInX0);
          // apply on the covariance matrix
          multiscattering(
              state, stepper, state.stepping, sigmaScat, mInteraction);
        }

        // Apply the Energy loss
        if (energyLoss) {
          // Get the material
          const Material& mat = mProperties.material();
          // Calculate gamma
          const double lgamma = E / m;
          // Energy loss and straggling - per unit length
          std::pair<double, double> eLoss = process_ionisationloss.dEds(
              m, lbeta, lgamma, mat, 1. * units::_mm);
          // Apply the energy loss
          const double dEdl   = state.stepping.navDir * eLoss.first;
          const double dE     = mProperties.thickness() * dEdl;
          const double sigmaE = state.stepping.navDir * eLoss.second;
          // Screen output
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Energy loss calculated to " << dE / 1_GeV << " GeV";
            return dstream.str();
          });
          energyloss(state,
                     stepper,
                     state.stepping,
                     p,
                     m,
                     E,
                     dE,
                     sigmaE,
                     mInteraction);
        }

        // This doesn't cost anything - do it regardless
        result.materialInX0 += mProperties.thicknessInX0();
        result.materialInL0 += mProperties.thicknessInL0();

        // Record the material interaction if configured to do so
        if (recordInteractions) {
          mInteraction.position = stepper.position(state.stepping);
          mInteraction.direction = stepper.direction(state.stepping);
          mInteraction.materialProperties = mProperties;
          mInteraction.pathCorrection = pCorrection;
          result.materialInteractions.push_back(std::move(mInteraction));
        }
      }
    }
  }

private:
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
  /// @param logAction is a callable function that returns a stremable object
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
