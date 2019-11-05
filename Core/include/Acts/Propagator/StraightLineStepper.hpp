// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Propagator/detail/StepperReturnState.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Propagator/StepperState.hpp"

namespace Acts {
  
/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
class StraightLineStepper {
 public:

  using Corrector = VoidIntersectionCorrector;
  using Jacobian = BoundMatrix;
  using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  using State = StepperState;

  /// Constructor
  StraightLineStepper() = default;

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D getField(State& /*state*/, const Vector3D& /*pos*/) const {
    // get the field from the cell
    return Vector3D(0., 0., 0.);
  }

  /// Global particle position accessor
  Vector3D position(const State& state) const { return state.pos; }

  /// Momentum direction accessor
  Vector3D direction(const State& state) const { return state.dir; }

  /// Momentum accessor
  double momentum(const State& state) const { return state.p; }

  /// Charge access
  double charge(const State& state) const { return state.q; }

  /// Time access
  double time(const State& state) const { return state.t0 + state.dt; }

  /// Tests if the state reached a surface
  ///
  /// @param [in] state State that is tests
  /// @param [in] surface Surface that is tested
  ///
  /// @return Boolean statement if surface is reached by state
  bool surfaceReached(const State& state, const Surface* surface) const {
    return surface->isOnSurface(state.geoContext, position(state),
                                direction(state), true);
  }

  /// @brief Final state builder without a target surface
  ///
  /// @tparam start_parameters_t Type of the start parameters
  /// @tparam end_parameters_t Type of the end parameters
  ///
  /// @param [in, out] state State of the propagation
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed,
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return std::tuple conatining the final state parameters, the jacobian & the accumulated path
  template<bool start_local, typename end_parameters_t>
  auto 
  buildState(State& state, bool reinitialize) const
  {	  
	  // The return type
	  using return_type = detail::return_state_type<start_local, end_parameters_t>;
	  // If the result should be local it is curvilinear
	  if constexpr (end_parameters_t::is_local_representation)
	  {
		 return covTransport.curvilinearState<return_type>(state, reinitialize);
	  }
	  // else it is free
	  else
	  {
		   return covTransport.freeState<return_type>(state, reinitialize);
	  }
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed,
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return std::tuple conatining the final state parameters, the jacobian & the accumulated path
  template<bool start_local>
  auto 
  buildState(State& state, const Surface& surface, bool reinitialize) const {
	  using return_type = detail::return_state_type<start_local, BoundParameters, Surface>;
	  	  return covTransport.boundState<return_type>(state, surface, reinitialize);
  }

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  void update(State& state, const BoundParameters& pars) const {
    const auto& mom = pars.momentum();
    state.pos = pars.position();
    state.dir = mom.normalized();
    state.p = mom.norm();
    state.dt = pars.time();

    if (pars.covariance().has_value()) {	
		  state.cov = *pars.covariance();
    }
  }

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3D& uposition,
              const Vector3D& udirection, double up, double time) const {
    state.pos = uposition;
    state.dir = udirection;
    state.p = up;
    state.dt = time;
  }

  /// Return a corrector
  VoidIntersectionCorrector corrector(State& /*state*/) const {
    return VoidIntersectionCorrector();
  }

  /// Perform a straight line propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///                The state contains the desired step size,
  ///                it can be negative during backwards track propagation,
  ///                and since we're using an adaptive algorithm, it can
  ///                be modified by the stepper class during propagation.
  ///
  /// @return the step size taken
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const {
    // use the adjusted step size
    const auto h = state.stepping.stepSize;
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = std::hypot(1., state.options.mass / state.stepping.p);
    // Update the track parameters according to the equations of motion
    state.stepping.pos += h * state.stepping.dir;
    state.stepping.dt += h * dtds;
    // Propagate the jacobian
    if (state.stepping.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = ActsSymMatrixD<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * state.options.mass * state.options.mass * state.stepping.q /
                (state.stepping.p * dtds);
      // Set the derivative factor the time
      state.stepping.derivative(3) = dtds;
      // Update jacobian and derivative
      state.stepping.jacTransport = D * state.stepping.jacTransport;
      state.stepping.derivative.template head<3>() = state.stepping.dir;
    }
    // state the path length
    state.stepping.pathAccumulated += h;

    // return h
    return h;
  }

private:
	/// The covariance transporter engine
	CovarianceTransport covTransport;
};
}  // namespace Acts
