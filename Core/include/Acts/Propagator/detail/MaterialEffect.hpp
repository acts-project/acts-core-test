// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/DefaultParameterDefinitions.hpp"

namespace Acts {

namespace detail {

  struct MultipleScattering
  {
    template <typename propagator_state_t,
              typename stepper_t,
              typename component_state_t,
              typename material_interaction_t>
    void
    operator()(propagator_state_t& /*no use*/,
               const stepper_t&        stepper,
               component_state_t&      component_state,
               double                  sigma,
               material_interaction_t& interaction) const
    {
      double sinTheta
          = std::sin(VectorHelpers::theta(stepper.direction(component_state)));
      double sigmaDeltaPhiSq   = sigma * sigma / (sinTheta * sinTheta);
      double sigmaDeltaThetaSq = sigma * sigma;
      // Record the material interaction
      interaction.sigmaPhi2   = sigmaDeltaPhiSq;
      interaction.sigmaTheta2 = sigmaDeltaThetaSq;
      // Good in any case for positive direction
      if (component_state.navDir == forward) {
        // Just add the multiple scattering component
        component_state.cov(ePHI, ePHI)
            += component_state.navDir * sigmaDeltaPhiSq;
        component_state.cov(eTHETA, eTHETA)
            += component_state.navDir * sigmaDeltaThetaSq;
      } else {
        // We check if the covariance stays positive
        double sEphi   = component_state.cov(ePHI, ePHI);
        double sEtheta = component_state.cov(eTHETA, eTHETA);
        if (sEphi > sigmaDeltaPhiSq && sEtheta > sigmaDeltaThetaSq) {
          // Noise removal is not applied if covariance would fall below 0
          component_state.cov(ePHI, ePHI) -= sigmaDeltaPhiSq;
          component_state.cov(eTHETA, eTHETA) -= sigmaDeltaThetaSq;
        }
      }
    }
  };

  ///
  struct EnergyLoss
  {
    template <typename propagator_state_t,
              typename stepper_t,
              typename component_state_t,
              typename material_interaction_t>
    void
    operator()(propagator_state_t& /*nouse*/,
               const stepper_t&        stepper,
               component_state_t&      component_state,
               const double            p,
               const double            m,
               const double            E,
               const double            dE,
               const double            deltaCov,
               material_interaction_t& interaction) const
    {
      // Check for energy conservation, and only apply momentum change
      // when kinematically allowed
      if (E + dE > m) {
        const double newP = std::sqrt((E + dE) * (E + dE) - m * m);
        component_state.p = std::copysign(newP, component_state.p);
        // Record the deltaP
        interaction.deltaP = p - newP;
        // Update the state/momentum
        stepper.update(component_state,
                       stepper.position(component_state),
                       stepper.direction(component_state),
                       std::copysign(newP, stepper.momentum(component_state)));
      }
      // Transfer this into energy loss straggling and apply to
      // covariance:
      // do that even if you had not applied energy loss due to
      // the kineamtic limit to catch the cases of deltE < MOP/MPV
      if (component_state.covTransport) {
        // current variance varity : 0
        const double sigmaQoverP = deltaCov / (p / E * p * p);
        // Save the material interaction
        interaction.sigmaQoP2 = sigmaQoverP * sigmaQoverP;
        // good in any case for positive direction
        if (component_state.navDir == forward) {
          component_state.cov(eQOP, eQOP)
              += component_state.navDir * sigmaQoverP * sigmaQoverP;
        } else {
          // check that covariance entry doesn't become negative
          double sEqop = component_state.cov(eQOP, eQOP);
          if (sEqop > sigmaQoverP * sigmaQoverP) {
            component_state.cov(eQOP, eQOP)
                += component_state.navDir * sigmaQoverP * sigmaQoverP;
          }
        }
      }
    }
  };
}
}
