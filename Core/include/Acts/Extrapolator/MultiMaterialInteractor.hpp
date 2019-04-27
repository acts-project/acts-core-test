// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/detail/EmptyEffects.hpp"
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

using InteractionPointVec = std::vector<Acts::MaterialInteraction>;

/// The Material interactor struct
///
/// This is a plugin to the Propagator that
/// performs material interaction on the currentSurface
/// of the Propagagor state
struct MultiMaterialInteractor
{
  /// Configuration for this MultiMaterialInteractor

  /// multiple scattering switch on/off
  bool multipleScattering = true;
  /// The scattering formula struct
  detail::HighlandScattering process_scattering;
  /// apply delta theta/phi on the covariance matrix
  detail::MultipleScattering multiplescattering;

  /// Empty Bethe-Heitler struct
  /// currently copy one component into multi equal components
  detail::EmptyEffect emptyEffect;
  /// apply dE and delta covariance on the matrix
  detail::EnergyLoss energyloss;

  /// Record material in detail
  bool recordInteractions = false;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    /// This one is only filled when recordInteractions is switched on
    /// record all the interactionPoints on surface
    std::map<const Surface*, InteractionPointVec> multiMaterialInteractions;
    /// the number of components
    int numComponents = 0;
  };

  using result_type = this_result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// only contains a split function currently
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// the energy loss part are assumed to be Bethe-Heitler function, which
  /// allows to split
  /// one component into several components with (weight,deltaE,variance),
  /// currently the Bethe-Heitler does nothing but split into several same
  /// components,
  /// to test if they collect same results in the multi-stepper. Then the
  /// multiple scattering part is considered.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  ///
  /// @to do Add Bethe-Heitler effect.
  ///
  template <typename propagator_state_t, typename stepper_t>
  void
  operator()(propagator_state_t& state,
             const stepper_t&    stepper,
             result_type&        result) const
  {
    debugLog(state, [&] { return std::string("in MultiMaterialInteractor."); });

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    // If switched off, then return - alows run-time configuration
    if (!multipleScattering && !recordInteractions) {
      return;
    }

    // A current surface has been already assigned by the navigator
    // check for material
    if (state.navigation.currentSurface
        && state.navigation.currentSurface->surfaceMaterial()) {
      // Let's set the pre/full/post update stage
      MaterialUpdateStage mStage = fullUpdate;
      // We are at the start surface
      if (state.navigation.startSurface == state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on start surface: post-update mode.");
        });
        mStage = postUpdate;
        // Or is it the target surface ?
      } else if (state.navigation.targetSurface
                 == state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on target surface: pre-update mode");
        });
        mStage = preUpdate;
      } else {
        debugLog(state, [&] {
          return std::string("Update while pass through: full mode.");
        });
      }

      /// Get the surface material & properties from them and continue if you
      /// found some
      const ISurfaceMaterial* sMaterial
          = state.navigation.currentSurface->surfaceMaterial();

      const double        m = state.options.mass;
      InteractionPointVec materialInteractionVec;

      /// typename of state column in the multi-state :
      /// std::list<tuple<state,weight,status>>
      using stateColType = decltype(state.stepping.stateCol);

      /// loop the single components in the list
      typename stateColType::iterator it = state.stepping.stateCol.begin();
      while (it != state.stepping.stateCol.end()) {
        // get the current single component state
        auto&        singlestate = std::get<0>(*it);
        double       weight      = std::get<1>(*it);
        const auto   status      = std::get<2>(*it);
        const double p           = singlestate.p;
        const double E           = std::sqrt(p * p + m * m);
        const double lbeta       = p / E;

        MaterialProperties mProperties = sMaterial->materialProperties(
            stepper.position(singlestate), singlestate.navDir, mStage);
        // Material properties (non-zero) have been found for this configuration
        if (mProperties) {
          /// more debugging output to the screen
          debugLog(state, [&] {
            return std::string("Material properties found for this surface.");
          });
          // Calculate the path correction
          double pCorrection = state.navigation.currentSurface->pathCorrection(
              state.geoContext,
              stepper.position(singlestate),
              stepper.direction(singlestate));

          // Scale the material properties
          mProperties *= pCorrection;
          const double tInX0 = mProperties.thicknessInX0();

          // Create the material interaction class, in case we record afterwards
          // Record the material interaction if configured to do so
          Acts::MaterialInteraction mInteraction;

          // To integrate process noise, we need to transport
          // the covariance to the current position in space
          // the 'true' indicates re-initializaiton of the further transport
          if (singlestate.covTransport) {
            stepper.covarianceTransport(singlestate, true);
          }

          // @brief if meets material surface , split the current single state
          // to get a new list of components.
          //
          // This should be according to the Bethe-Heitler pdf
          // in the function, return the list of multicomponents, each component
          // carries (weight,dmean,dvariance)
          // weight represents the weight of each newly created component
          // dmean and dvariance represents the delta Energyloss and the delta
          // variance
          // from the newly created component.
          //
          // @note current not use the Bethe-Heitler pdf, just make a list of
          // copied component, to see if they act equally in the multi-stepper,
          // the number of split is set to 2

          // the mixture represents the vector of (weight,dmean,dvariance)
          // struct
          auto mixture = emptyEffect.getMixture(tInX0, p);
          for (const auto& mix : mixture) {
            // energy loss for each created component
            const double dE        = mix.mean;
            const double sigmaE    = mix.variance;
            const double pdfWeight = mix.weight;
            energyloss(
                state, stepper, singlestate, p, m, E, dE, sigmaE, mInteraction);

            // multiple scattering
            // currently multiple scattering is not considered a Bethe-Heitler
            // process
            if (multipleScattering && singlestate.covTransport) {
              double sigmaScat = process_scattering(p, lbeta, tInX0);
              multiplescattering(
                  state, stepper, singlestate, sigmaScat, mInteraction);
            }

            // get the weight of the new component
            double newWeight = pdfWeight * weight;
            state.stepping.stateCol.push_front(
                std::make_tuple(std::move(singlestate), newWeight, status));
          }

          // delete the current component
          it = state.stepping.stateCol.erase(it);

          if (recordInteractions) {
            mInteraction.surface            = state.navigation.currentSurface;
            mInteraction.position           = stepper.position(singlestate);
            mInteraction.direction          = stepper.direction(singlestate);
            mInteraction.materialProperties = mProperties;
            mInteraction.pathCorrection     = pCorrection;
            materialInteractionVec.push_back(std::move(mInteraction));
          }  // end of record
        }
      }
      // record the material interaction
      if (recordInteractions) {
        result.multiMaterialInteractions.insert(
            std::pair<const Surface*, InteractionPointVec>(
                state.navigation.currentSurface,
                std::move(materialInteractionVec)));
      }
      // record the number of components at the last step
      result.numComponents = state.stepping.stateCol.size();
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
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
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

}  // end of namespace Acts
