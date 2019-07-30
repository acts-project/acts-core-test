// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename B, typename C, typename E, typename A>
Acts::Vector3D Acts::MultiEigenStepper<B, C, E, A>::position(
    const State& state) const {
  Vector3D pos = Vector3D(0, 0, 0);
  for (const auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    Vector3D component_pos =
        EigenStepperType::position(std::get<0>(tuple_state));
    pos += component_pos * std::get<1>(tuple_state);
  }
  return pos;
}

template <typename B, typename C, typename E, typename A>
Acts::Vector3D Acts::MultiEigenStepper<B, C, E, A>::direction(
    const State& state) const {
  Vector3D dir = Vector3D(0, 0, 0);
  for (const auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    Vector3D component_dir =
        EigenStepperType::direction(std::get<0>(tuple_state));
    dir += component_dir * std::get<1>(tuple_state);
  }
  return dir;
}

template <typename B, typename C, typename E, typename A>
double Acts::MultiEigenStepper<B, C, E, A>::momentum(const State& state) const {
  double mom = 0.;
  for (const auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    double component_mom = EigenStepperType::momentum(std::get<0>(tuple_state));
    mom += component_mom * std::get<1>(tuple_state);
  }
  return mom;
}

template <typename B, typename C, typename E, typename A>
double Acts::MultiEigenStepper<B, C, E, A>::time(const State& state) const {
  double total_time = 0.;
  for (const auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    double component_time = EigenStepperType::time(std::get<0>(tuple_state));
    total_time += component_time * std::get<1>(tuple_state);
  }
  return total_time;
}

template <typename B, typename C, typename E, typename A>
bool Acts::MultiEigenStepper<B, C, E, A>::surfaceReached(
    State& state, const Surface* surface) const {
  // status is true when there all free components are on surface
  bool status = true;
  for (auto& tuple_state : state.stateCol) {
    auto& singlestate = std::get<0>(tuple_state);
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    // if this component is not on surface, status is set false
    if (!surface->isOnSurface(singlestate.geoContext,
                              EigenStepperType::position(singlestate),
                              EigenStepperType::direction(singlestate), true)) {
      status = false;
    }
    // if this component on surface: locked
    else {
      std::get<2>(tuple_state) = StateStatus::Locked;
    }
  }
  if (status == false) {
    return false;
  } else {
    // if all the components are on the surface(except the dead ones), set them
    // free and set the currentSurface in Navigator.
    releaseComponents(state);
    return true;
  }
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::releaseComponents(
    State& state) const {
  for (auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) == StateStatus::Locked) {
      std::get<2>(tuple_state) = StateStatus::Free;
    }
  }
}

template <typename B, typename C, typename E, typename A>
template <typename options_t>
auto Acts::MultiEigenStepper<B, C, E, A>::targetSurface(
    State& state, const Surface* surface, const options_t& navOpts,
    const Corrector& navCorr) const -> SurfaceIntersection {
  SurfaceIntersection minimumIntersection;
  double minDist = std::numeric_limits<double>::max();
  for (auto& tuple_state : state.stateCol) {
    /// only target the free component
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    auto& singlestate = std::get<0>(tuple_state);
    auto target =
        EigenStepperType::targetSurface(singlestate, surface, navOpts, navCorr);
    if (target && minDist > target.intersection.pathLength) {
      minDist = target.intersection.pathLength;
      minimumIntersection = target;
    }
    /// a protection to avoid the flip components,
    /// which is an abnormal propagate
    /// the dead case should be considered in the Fitter best,
    /// currently we simply set it dead here
    if (direction(state).dot(direction(singlestate)) < 0 || !target) {
      std::get<2>(tuple_state) = StateStatus::Dead;
    }
  }
  /// delete the dead component
  deleteComponents(state);
  return minimumIntersection;
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::normalizeComponents(
    State& state) const {
  double weight_sum = 0;
  for (const auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) == StateStatus::Dead)
      continue;
    weight_sum += std::get<1>(tuple_state);
  }
  for (auto& tuple_state : state.stateCol) {
    if (std::get<2>(tuple_state) == StateStatus::Dead)
      continue;
    std::get<1>(tuple_state) = std::get<1>(tuple_state) / weight_sum;
  }
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::deleteComponents(State& state) const {
  auto& col = state.stateCol;
  typename std::list<std::tuple<SingleStateType, double, StateStatus>>::iterator
      it = col.begin();
  while (it != col.end()) {
    if (std::get<2>(*it) == StateStatus::Dead) {
      it = col.erase(it);
    } else {
      ++it;
    }
  }
  /// normalize the weight to 1.
  normalizeComponents(state);
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::updateStepSize(
    State& state, const Corrector& navCorr, double navigationStep,
    bool release) const {
  for (auto& tuple_state : state.stateCol) {
    // only deal with the free
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    auto& singlestate = std::get<0>(tuple_state);
    EigenStepperType::updateStepSize(singlestate, navCorr, navigationStep,
                                     release);
  }
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::releaseStep(State& state,
                                                      cstep::Type type) const {
  for (auto& tuple_state : state.stateCol) {
    // only deal with the free
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    auto& singlestate = std::get<0>(tuple_state);
    EigenStepperType::releaseStep(singlestate, type);
  }
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::updateStepSize(
    State& state, double abortStep, cstep::Type type) const {
  for (auto& tuple_state : state.stateCol) {
    // only deal with the free
    if (std::get<2>(tuple_state) != StateStatus::Free)
      continue;
    auto& singlestate = std::get<0>(tuple_state);
    EigenStepperType::updateStepSize(singlestate, abortStep, type);
  }
  state.stepSize.update(abortStep, cstep::aborter);
}

template <typename B, typename C, typename E, typename A>
auto Acts::MultiEigenStepper<B, C, E, A>::boundState(State& state,
                                                     const Surface& surface,
                                                     bool /*unused*/) const
    -> BoundState {
  // covariance no use in the mcs
  std::unique_ptr<const Covariance> covPtr = nullptr;
  // Create the bound parameters
  BoundParameters parameters((std::get<0>(*state.stateCol.begin())).geoContext,
                             std::move(covPtr), state.pos, state.p * state.dir,
                             state.q,
                             (std::get<0>(*state.stateCol.begin())).t0 +
                                 (std::get<0>(*state.stateCol.begin())).dt,
                             surface.getSharedPtr());
  // Create the bound state
  BoundState bState{std::move(parameters), Jacobian::Identity(),
                    state.pathAccumulated};
  /// Return the State
  return bState;
}

template <typename B, typename C, typename E, typename A>
auto Acts::MultiEigenStepper<B, C, E, A>::curvilinearState(
    State& state, bool /*unused*/) const -> CurvilinearState {
  // covariance no use in the mcs
  std::unique_ptr<const Covariance> covPtr = nullptr;
  // Create the curvilinear parameters
  CurvilinearParameters parameters(std::move(covPtr), state.pos,
                                   state.p * state.dir, state.q,
                                   std::get<0>(*state.stateCol.begin()).t0 +
                                       std::get<0>(*state.stateCol.begin()).dt);
  // Create the bound state
  CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                             state.pathAccumulated};
  /// Return the State
  return curvState;
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::update(
    SingleStateType& singlestate, const BoundParameters& pars) const {
  EigenStepperType::update(singlestate, pars);
}

template <typename B, typename C, typename E, typename A>
void Acts::MultiEigenStepper<B, C, E, A>::update(SingleStateType& singlestate,
                                                 const Vector3D& uposition,
                                                 const Vector3D& udirection,
                                                 double up,
                                                 double component_time) const {
  EigenStepperType::update(singlestate, uposition, udirection, up,
                           component_time);
}

template <typename B, typename C, typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double> Acts::MultiEigenStepper<B, C, E, A>::step(
    propagator_state_t& state) const {
  // a variable to record the combination of the pathlength of the compact
  double combinedPath = 0;
  // loop all the components that are free
  for (auto& tuple_state : state.stepping.stateCol) {
    /// if the status is locked or dead, ignore it
    if (std::get<2>(tuple_state) != StateStatus::Free) {
      continue;
    }
    auto& singlestate = std::get<0>(tuple_state);

    // Runge-Kutta integrator state
    auto& sd = singlestate.stepData;

    double h2, half_h;
    double error_estimate;

    // First Runge-Kutta point (at current position)
    sd.B_first = EigenStepperType::getField(singlestate, singlestate.pos);
    if (!singlestate.extension.validExtensionForStep(state, *this,
                                                     singlestate) ||
        !singlestate.extension.k1(state, *this, singlestate, sd.k1,
                                  sd.B_first)) {
      return 0.;
    }

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the
    // local
    // integration error. The results are stated in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory
    const auto tryRungeKuttaStep = [&](const double h) -> bool {
      // State the square and half of the step size
      h2 = h * h;
      half_h = h * 0.5;

      // Second Runge-Kutta point
      const Vector3D pos1 =
          singlestate.pos + half_h * singlestate.dir + h2 * 0.125 * sd.k1;
      sd.B_middle = EigenStepperType::getField(singlestate, pos1);
      if (!singlestate.extension.k2(state, *this, singlestate, sd.k2,
                                    sd.B_middle, half_h, sd.k1)) {
        return false;
      }

      // Third Runge-Kutta point
      if (!singlestate.extension.k3(state, *this, singlestate, sd.k3,
                                    sd.B_middle, half_h, sd.k2)) {
        return false;
      }

      // Last Runge-Kutta point
      const Vector3D pos2 =
          singlestate.pos + h * singlestate.dir + h2 * 0.5 * sd.k3;
      sd.B_last = EigenStepperType::getField(singlestate, pos2);
      if (!singlestate.extension.k4(state, *this, singlestate, sd.k4, sd.B_last,
                                    h, sd.k3)) {
        return false;
      }

      // Return an estimate of the local integration error
      error_estimate = std::max(
          h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>(), 1e-20);
      return (error_estimate <= state.options.tolerance);
    };

    double stepSizeScaling;

    // Select and adjust the appropriate Runge-Kutta step size as given
    // ATL-SOFT-PUB-2009-001
    while (!tryRungeKuttaStep(singlestate.stepSize)) {
      stepSizeScaling =
          std::min(std::max(0.25, std::pow((state.options.tolerance /
                                            std::abs(2. * error_estimate)),
                                           0.25)),
                   4.);
      if (stepSizeScaling == 1.) {
        break;
      }
      singlestate.stepSize = singlestate.stepSize * stepSizeScaling;

      // If step size becomes too small the particle remains at the initial
      // place
      if (singlestate.stepSize * singlestate.stepSize <
          state.options.stepSizeCutOff * state.options.stepSizeCutOff) {
        // Not moving due to too low momentum needs an aborter
        return EigenStepperError::StepSizeStalled;
      }
    }

    // use the adjusted step size
    const double h = singlestate.stepSize;

    // When doing error propagation, update the associated Jacobian matrix
    if (singlestate.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D;
      if (!singlestate.extension.finalize(state, *this, singlestate, h, D)) {
        return EigenStepperError::StepInvalid;
      }

      // for moment, only update the transport part
      singlestate.jacTransport = D * singlestate.jacTransport;
    } else {
      if (!singlestate.extension.finalize(state, *this, singlestate, h)) {
        return EigenStepperError::StepInvalid;
      }
    }

    // Update the track parameters according to the equations of motion
    singlestate.pos += h * singlestate.dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
    singlestate.dir += h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
    singlestate.dir /= singlestate.dir.norm();
    singlestate.derivative.template head<3>() = singlestate.dir;
    singlestate.derivative.template segment<3>(3) = sd.k4;
    singlestate.pathAccumulated += h;

    // for the multistepper simply record the pathlength with combination of
    // components
    //
    state.stepping.pathAccumulated += h * std::get<1>(tuple_state);
    combinedPath += h * std::get<1>(tuple_state);
  }
  // the pos/dir/mom after RKN is a combination of the all
  state.stepping.pos = position(state.stepping);
  state.stepping.dir = direction(state.stepping);
  state.stepping.p = momentum(state.stepping);
  return combinedPath;
}
