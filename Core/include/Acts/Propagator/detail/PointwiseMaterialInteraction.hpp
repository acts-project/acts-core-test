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
 
 /// @brief This function evaluates the material properties
 /// 
 /// @tparam propagator_state_t Type of the propagator state
 /// @tparam stepper_t Type of the stepper
 /// @param [in, out] state State of the propagation
 /// @param [in] stepper Stepper of the propagation
 ///
 /// @return Material properties for the interaction
 template <typename propagator_state_t, typename stepper_t>
 MaterialProperties evaluateMaterialProperties(propagator_state_t& state, const stepper_t& stepper)
 {
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
	  return sMaterial->materialProperties(
		  stepper.position(state.stepping), state.stepping.navDir, mStage);
   }
 
  /// @brief This function evaluates the multiple scattering
 /// 
 /// @tparam propagator_state_t Type of the propagator state
 /// @tparam stepper_t Type of the stepper
 /// @param [in, out] state State of the propagation
 /// @param [in] stepper Stepper of the propagation
 /// @param [in] tInX0 Thickness in X0 of the material
 /// @param [in] p Momentum of the particle
 /// @param [in] lbeta Beta value of the particle
 ///
 /// @return Pair containing the modification of the variance in phi and theta due to multiple scattering
  template <typename propagator_state_t, typename stepper_t>
  std::pair<double, double> evaluateMultipleScattering(propagator_state_t& state, const stepper_t& stepper, double tInX0, double p, double lbeta)
  {
	  // Retrieve the scattering contribution
	  bool isElectron = state.options.absPdgCode == 11;
	  const double sigmaScat = detail::HighlandScattering::sigmaAngle(p, lbeta, tInX0, isElectron);
	  const double sinTheta =
		  std::sin(VectorHelpers::theta(stepper.direction(state.stepping)));
	  const double sigmaDeltaPhiSq =
		  sigmaScat * sigmaScat / (sinTheta * sinTheta);
	  const double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
	  // Good in any case for positive direction
	  if (state.stepping.navDir == forward) {
		// Just add the multiple scattering component
		state.stepping.cov(ePHI, ePHI) +=
			state.stepping.navDir * sigmaDeltaPhiSq;
		state.stepping.cov(eTHETA, eTHETA) +=
			state.stepping.navDir * sigmaDeltaThetaSq;
	  } else {
		// We check if the covariance stays positive
		const double sEphi = state.stepping.cov(ePHI, ePHI);
		const double sEtheta = state.stepping.cov(eTHETA, eTHETA);
		if (sEphi > sigmaDeltaPhiSq && sEtheta > sigmaDeltaThetaSq) {
		  // Noise removal is not applied if covariance would fall below 0
		  state.stepping.cov(ePHI, ePHI) -= sigmaDeltaPhiSq;
		  state.stepping.cov(eTHETA, eTHETA) -= sigmaDeltaThetaSq;
		}
	  }
	  return std::make_pair(sigmaDeltaPhiSq, sigmaDeltaThetaSq);
  }
   
 /// @brief This function evaluates the covariance contribution due to energy loss in material
 /// 
 /// @tparam propagator_state_t Type of the propagator state
 /// @param [in, out] state State of the propagation
 /// @param [in] thickness Thickness of the material
 /// @param [in] eLoss Width of the energy loss in the material
 /// @param [in] p Momentum of the particle
 /// @param [in] lbeta Beta value of the particle
 ///
 /// @return Contribution of the energy loss to the covariance matrix
  template <typename propagator_state_t>
  double evaluateEnergyLossForCovariance(propagator_state_t& state, double thickness, double eLoss, double p, double lbeta)
  {
	// Calculate the straggling
	const double sigmaQoverP =
		thickness * eLoss / (lbeta * p * p);
	// Save the material interaction
	const double sigmaQoverP2 = sigmaQoverP * sigmaQoverP;
	// Good in any case for positive direction
	if (state.stepping.navDir == forward) {
	  state.stepping.cov(eQOP, eQOP) +=
		  state.stepping.navDir * sigmaQoverP2;
	} else {
	  // Check that covariance entry doesn't become negative
	  const double sEqop = state.stepping.cov(eQOP, eQOP);
	  if (sEqop > sigmaQoverP2) {
		state.stepping.cov(eQOP, eQOP) +=
			state.stepping.navDir * sigmaQoverP2;
	  }
	}
	return sigmaQoverP2;
  }
  
 /// @brief This function evaluates the effect on a particle due to energy loss in material
 /// 
 /// @tparam propagator_state_t Type of the propagator state 
 /// @tparam stepper_t Type of the stepper
 /// @param [in, out] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
 /// @param [in] mProperties Properties of the material
 /// @param [in] E Energy of the particle
 /// @param [in] p Momentum of the particle
 /// @param [in] m Mass of the particle
 /// @param [in] lbeta Beta value of the particle
 ///
 /// @return Pair containing the momentum change and the covariance contribution due to the energy loss
   template <typename propagator_state_t, typename stepper_t>
  std::pair<double, double> evaluateEnergyLoss(propagator_state_t& state, const stepper_t& stepper, const MaterialProperties& mProperties, double E, double p, double m, double lbeta)
  {
	 using namespace Acts::UnitLiterals;
	 
     // Get the material
	  const Material& mat = mProperties.material();
	  // Calculate gamma
	  const double lgamma = E / m;
	  // Energy loss and straggling - per unit length
	  std::pair<double, double> eLoss =
		  detail::IonisationLoss::dEds(m, lbeta, lgamma, mat, 1_mm);
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
		// Update the state/momentum
		stepper.update(
			state.stepping, stepper.position(state.stepping),
			stepper.direction(state.stepping),
			std::copysign(newP, stepper.momentum(state.stepping)),
			stepper.time(state.stepping));
	  }
	  // Transfer this into energy loss straggling and apply to
	  // covariance:
	  // do that even if you had not applied energy loss due to
	  // the kineamtic limit to catch the cases of deltE < MOP/MPV
	  if (state.stepping.covTransport) {
		  result.second = evaluateEnergyLossForCovariance(state, mProperties.thickness(), eLoss.second, p, lbeta);
	  }
	  return result;
  }
  
   /// @brief This function evaluates the effect of material on a particle
 /// 
 /// @tparam propagator_state_t Type of the propagator state 
 /// @tparam stepper_t Type of the stepper
 /// @param [in, out] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
 /// @param [in] mProperties Properties of the material
 ///
 /// @return Array containing the changes to the variance of phi, theta, the momentum and the variance of QoP
  template <typename propagator_state_t, typename stepper_t>
  std::array<double, 4> evaluateMaterialInteraction(propagator_state_t& state, const stepper_t& stepper, const MaterialProperties& mProperties, bool multipleScattering = true, bool energyLoss = true)
  {
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
		  sigmaAngles = evaluateMultipleScattering(state, stepper, mProperties.thicknessInX0(), p, lbeta);
        }

        // Apply the Energy loss
        if (energyLoss) {
			sigmaQoP = evaluateEnergyLoss(state, stepper, mProperties, E, p, m, lbeta);
        }
        return {sigmaAngles.first, sigmaAngles.second, sigmaQoP.first, sigmaQoP.second};
  }
}
}  // end of namespace Acts