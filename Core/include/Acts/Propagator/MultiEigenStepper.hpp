// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/EigenStepper.hpp"

namespace Acts {

enum class StateStatus { Free = 0, Locked = 1, Dead = 2 };

/// @brief multicomponent stepper(mcs) is based on the Single Stepper
/// implementation developed for Gaussian Sum Filter (GSF)
///
/// the state of MC contains a list of single state
/// each component contains its own pos, dirand stepSize
/// in the step() method,
/// loop all the single components and do caculating as in
/// the single stepper.
/// In surfaceReached() method, determine if all components
/// are on the current surface.
/// In targetSurface() method, collect the candidate surfaces in Navigator
/// with the combination of components (pos,dir),
/// and update the stepSize of each components
///
/// in mcs, each single component owns a status:
/// Free - not on surface
/// Lock - on surface
/// Dead - can not target the surface
///
/// @to do move the job of deleting components to the Fitting stage
/// @to do compliant with the MultiParameter class

template <typename BField, typename corrector_t = VoidIntersectionCorrector,
          typename extensionlist_t = StepperExtensionList<DefaultExtension>,
          typename auctioneer_t = detail::VoidAuctioneer>
class MultiEigenStepper
    : public EigenStepper<BField, corrector_t, extensionlist_t, auctioneer_t> {
 public:
  using cstep = detail::ConstrainedStep;
  using Corrector = corrector_t;

  /// Jacobian, Covariance and State defintions
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;

  /// @note Currently the BoundState/CurvilinearState is defined for Single
  /// component which is a combination of all components, the Jacobian for this
  /// meaning is nonsense
  /// This should be replaced by
  /// std::tuple<MultiBoundParameters,std::list<Jacobian>,double> or other wise
  /// structure with Jacobians
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;
  using EigenStepperType =
      EigenStepper<BField, corrector_t, extensionlist_t, auctioneer_t>;
  using SingleStateType = typename EigenStepperType::State;

  using SurfaceIntersection = ObjectIntersection<Surface>;

  /// @brief State for track parameter propagation
  ///
  /// the State behave like the SingleState in Navigator, while contains
  /// information of each single components
  struct State {
    /// Constructor from the initial track parameters
    /// construct the multi components from one single component
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direction w.r.t momentum
    /// @param [in] ssize is the maximum step size
    ///
    template <typename parameters_t>
    explicit State(std::reference_wrapper<const GeometryContext> gctx,
                   std::reference_wrapper<const MagneticFieldContext> mctx,
                   const parameters_t& par, NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
        : pos(par.position()),
          dir(par.momentum().normalized()),
          p(par.momentum().norm()),
          q(par.charge()),
          navDir(ndir),
          stepSize(ndir * std::abs(ssize)) {
      /// initialize the MC state with one component
      stateCol.push_back(
          std::make_tuple(SingleStateType(gctx, mctx, par, ndir, ssize), 1.,
                          StateStatus::Free));
      /// remember the start parameters
      startPos = pos;
      startDir = dir;
    }
    /// the list of <singleState, weight, status>
    std::list<std::tuple<SingleStateType, double, StateStatus>> stateCol;

    /// Global start particle position
    Vector3D startPos = Vector3D(0., 0., 0.);

    /// Momentum start direction (normalized)
    Vector3D startDir = Vector3D(1., 0., 0.);

    /// Global particle position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1., 0., 0.);

    /// Momentum
    double p = 0.;

    /// The charge
    /// Currently suppose every component has the same charge
    double q = 1.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// Currently a condition if all components should transport the covariance
    bool covTransport = false;

    /// Accummulated path length state with the combination of the alive
    /// components
    double pathAccumulated = 0.;

    /// Currently used in setting the navigation option
    cstep stepSize{std::numeric_limits<double>::max()};
  };

  /// @brief Global particle position accessor
  /// get the combination of the position of all Free components
  Vector3D position(const State& state) const;

  /// @brief Global direction accessor
  /// get the combination of the direction of all Free components
  Vector3D direction(const State& state) const;

  /// @brief momentum accessor
  /// @brief Global get the combination of the momentum of all Free components
  double momentum(const State& state) const;

  using EigenStepperType::direction;
  using EigenStepperType::getField;
  using EigenStepperType::momentum;
  using EigenStepperType::position;
  using EigenStepperType::time;

  /// Charge access
  double charge(const State& state) const { return state.q; }

  /// Time access
  double time(const State& state) const;

  /// Constructor requires knowledge of the detector's magnetic field
  MultiEigenStepper(BField bField = BField()) : EigenStepperType(bField) {}

  /// similar method of the Single stepper
  /// @brief getField from 1st single-component in the list
  Vector3D getField(State& state, const Vector3D& pos) const {
    return EigenStepperType::getField(std::get<0>(*state.stateCol.begin()),
                                      pos);
  }

  /// @brief Tests if all the single states reach a surface
  /// if all single states reach(except the dead ones) successfully,
  /// return true; otherwise return false;
  ///
  /// the single state that is successfully reached is set to locked
  /// then if all single states locked, free all of them
  ///
  /// @param [in] state State is the mcs state
  /// @param [in] surface Surface that is tested
  ///
  /// @return Boolean statement if surface is reached by state
  bool surfaceReached(State& state, const Surface* surface) const;

  /// @brief Release all the components that are Locked
  void releaseComponents(State& state) const;

  /// @brief this update the stepSize of all the components
  /// to the candidate surfaces/layers/boundaries in the Navigator
  /// @param [in] state State is the mcs
  /// @param [in] surface Surface that is tested
  /// @param [in] navigator options
  /// @param [in] navigator corrections
  ///
  /// @return Intersection if there is an intersection
  //
  /// @to do: the dead component should be taken into consideration in the
  /// Fitter stage
  template <typename options_t>
  SurfaceIntersection targetSurface(State& state, const Surface* surface,
                                    const options_t& navOpts,
                                    const Corrector& navCorr) const;

  /// @brief reweight the free components
  /// the free and locked components are reweighted to 1
  void normalizeComponents(State& state) const;

  /// @brief The method to delete the dead components
  /// this would be done at Fitting step - in the future
  void deleteComponents(State& state) const;

  /// @brief get a sinlge parameter of combination of multi component on a
  /// surface, the jocobian is nonsence here
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter would be MultiBoundParameters
  BoundState boundState(State& state, const Surface& surface,
                        bool reinitialize = true) const;

  /// @brief get a sinlge parameter of combination of multi component on a
  /// surface, the jocobian is nonsence here
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter would be MultiCurvilinearParameters
  CurvilinearState curvilinearState(State& state,
                                    bool reinitialize = true) const;

  /// updateStepSize method call at navigator
  /// use to udpate each stepsize with the combination value
  void updateStepSize(State& state, const Corrector& navCorr,
                      double navigationStep, bool release = false) const;

  void releaseStep(State& state, cstep::Type type = cstep::actor) const;

  /// call in StandardAborter,
  /// set a pathlimit of the combination component
  void updateStepSize(State& state, double abortStep,
                      cstep::Type type = cstep::aborter) const;

  /// @brief update for the single state, update singlestate to some parameters
  void update(SingleStateType& singlestate, const BoundParameters& pars) const;

  /// @brief update for the single state, update singlestate direction and p
  void update(SingleStateType& singlestate, const Vector3D& uposition,
              const Vector3D& udirection, double up, double time) const;

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const;
};
}  // namespace Acts

#include "Acts/Propagator/MultiEigenStepper.ipp"
