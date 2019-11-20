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
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/detail/InteractionFormulas.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace detail {
struct PointwiseMaterialInteraction {
  /// Default constructor
  PointwiseMaterialInteraction() = default;

  /// @brief This function evaluates the material properties
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in, out] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  ///
  /// @return Material properties for the interaction
  template <typename propagator_state_t, typename stepper_t>
  MaterialProperties evaluateMaterialProperties(
      const propagator_state_t& state, const stepper_t& stepper) const {
    // Let's set the pre/full/post update stage
    MaterialUpdateStage mStage = fullUpdate;
    // We are at the start surface
    if (state.navigation.startSurface == state.navigation.currentSurface) {
      mStage = postUpdate;
      // Or is it the target surface ?
    } else if (state.navigation.targetSurface ==
               state.navigation.currentSurface) {
      mStage = preUpdate;
    }

    // Get the surface material & properties from them and continue if you
    // found some
    const ISurfaceMaterial* sMaterial =
        state.navigation.currentSurface->surfaceMaterial();
    return sMaterial->materialProperties(stepper.position(state.stepping),
                                         state.stepping.navDir, mStage);
  }

  /// @brief This function evaluates the effect of material on a particle
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in, out] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] mProperties Properties of the material
  ///
  /// @return Array containing the changes to the variance of phi, theta, the
  /// momentum and the variance of QoP
  template <typename propagator_state_t, typename stepper_t>
  std::array<double, 4> evaluateMaterialInteraction(
      const propagator_state_t& state, const stepper_t& stepper,
      const MaterialProperties& mProperties, bool multipleScattering = true,
      bool energyLoss = true) const {
    // The momentum at current position
    const double p = stepper.momentum(state.stepping);
    const double m = state.options.mass;
    const double E = std::sqrt(p * p + m * m);
    const double lbeta = p / E;
    // Modifications due to material interactions
    std::pair<double, double> sigmaAngles;
    std::pair<double, double> sigmaQoP;

    // Apply the multiple scattering
    // - only when you do covariance transport
    if (multipleScattering && state.stepping.covTransport) {
      sigmaAngles = evaluateMultipleScattering(
          state, stepper, mProperties.thicknessInX0(), p, lbeta);
    }

    // Apply the Energy loss
    if (energyLoss) {
      sigmaQoP = evaluateEnergyLoss(state, mProperties, E, p, m, lbeta);
    }
    return {sigmaAngles.first, sigmaAngles.second, sigmaQoP.first,
            sigmaQoP.second};
  }

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

 private:
  /// @brief This function evaluates the multiple scattering
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] tInX0 Thickness in X0 of the material
  /// @param [in] p Momentum of the particle
  /// @param [in] lbeta Beta value of the particle
  ///
  /// @return Pair containing the modification of the variance in phi and theta
  /// due to multiple scattering
  template <typename propagator_state_t, typename stepper_t>
  std::pair<double, double> evaluateMultipleScattering(
      const propagator_state_t& state, const stepper_t& stepper, double tInX0,
      double p, double lbeta) const {
    // Retrieve the scattering contribution
    bool isElectron = state.options.absPdgCode == 11;
    const double sigmaScat = scattering(p, lbeta, tInX0, isElectron);
    const double sinTheta =
        std::sin(VectorHelpers::theta(stepper.direction(state.stepping)));
    const double sigmaDeltaPhiSq =
        sigmaScat * sigmaScat / (sinTheta * sinTheta);
    const double sigmaDeltaThetaSq = sigmaScat * sigmaScat;

    return std::make_pair(sigmaDeltaPhiSq, sigmaDeltaThetaSq);
  }

  /// @brief This function evaluates the covariance contribution due to energy
  /// loss in material
  ///
  /// @param [in] thickness Thickness of the material
  /// @param [in] eLoss Width of the energy loss in the material
  /// @param [in] p Momentum of the particle
  /// @param [in] lbeta Beta value of the particle
  ///
  /// @return Contribution of the energy loss to the covariance matrix
  double evaluateEnergyLossForCovariance(double thickness, double eLoss,
                                         double p, double lbeta) const {
    // Calculate the straggling
    const double sigmaQoverP = thickness * eLoss / (lbeta * p * p);
    // Return the material interaction
    return sigmaQoverP * sigmaQoverP;
  }

  /// @brief This function evaluates the effect on a particle due to energy loss
  /// in material
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @param [in] state State of the propagation
  /// @param [in] mProperties Properties of the material
  /// @param [in] E Energy of the particle
  /// @param [in] p Momentum of the particle
  /// @param [in] m Mass of the particle
  /// @param [in] lbeta Beta value of the particle
  ///
  /// @return Pair containing the momentum change and the covariance
  /// contribution due to the energy loss
  template <typename propagator_state_t>
  std::pair<double, double> evaluateEnergyLoss(
      const propagator_state_t& state, const MaterialProperties& mProperties,
      double E, double p, double m, double lbeta) const {
    using namespace Acts::UnitLiterals;

    // Get the material
    const Material& mat = mProperties.material();
    // Calculate gamma
    const double lgamma = E / m;
    // Energy loss and straggling - per unit length
    std::pair<double, double> eLoss =
        ionisationloss.dEds(m, lbeta, lgamma, mat, 1_mm);
    // Apply the energy loss
    const double dEdl = state.stepping.navDir * eLoss.first;
    const double dE = mProperties.thickness() * dEdl;
    std::pair<double, double> result;
    // Check for energy conservation, and only apply momentum change
    // when kinematically allowed
    if (E + dE > m) {
      // Calcuate the new momentum
      const double newP = std::sqrt((E + dE) * (E + dE) - m * m);
      // Record the deltaP
      result.first = p - newP;
    } else {
      result.first = 0.;
    }
    // Transfer this into energy loss straggling and apply to
    // covariance:
    // do that even if you had not applied energy loss due to
    // the kineamtic limit to catch the cases of deltE < MOP/MPV
    if (state.stepping.covTransport) {
      result.second = evaluateEnergyLossForCovariance(mProperties.thickness(),
                                                      eLoss.second, p, lbeta);
    }
    return result;
  }

  /// The scattering formula struct
  detail::HighlandScattering scattering;
  /// The energy loss formula struct
  detail::IonisationLoss ionisationloss;
};
}  // namespace detail
}  // end of namespace Acts