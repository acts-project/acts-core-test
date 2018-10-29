// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
<<<<<<< HEAD
#include "Acts/Propagator/StepperExtensionList.hpp"
=======
#include "Acts/Propagator/ExtensionList.hpp"
>>>>>>> Auctioneer setting in EigenStepper added; Covariance test for test cases a&b in unit test added
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @brief Runge-Kutta-Nystroem stepper based on Eigen implementation
/// for the following ODE:
///
/// r = (x,y,z)    ... global position
/// T = (Ax,Ay,Az) ... momentum direction (normalized)
///
/// dr/ds = T
/// dT/ds = q/p * (T x B)
///
/// with s being the arc length of the track, q the charge of the particle,
/// p its momentum and B the magnetic field
///
template <typename BField,
          typename corrector_t     = VoidCorrector,
          typename extensionlist_t = StepperExtensionList<DefaultExtension>,
          typename auctioneer_t    = detail::VoidAuctioneer>
class EigenStepper
{

private:
  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    using type = BoundParameters;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    using type = CurvilinearParameters;
  };

public:
  using cstep = detail::ConstrainedStep;

  /// @brief Storage of magnetic field and the sub steps during a RKN4 step
  struct StepData
  {
    /// Magnetic fields
    Vector3D B_first, B_middle, B_last;
    /// k_i of the RKN4 algorithm
    Vector3D k1, k2, k3, k4;
  };

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State
  {
    /// Constructor from the initial track parameters
    /// @param[in] par The track parameters at start
    /// @param[in] ndir The navigation direciton w.r.t momentum
    /// @param[in] sszice is the maximum step size
    ///
    /// @note the covariance matrix is copied when needed
    template <typename T>
    explicit State(const T&            par,
                   NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
      : pos(par.position())
      , dir(par.momentum().normalized())
      , p(par.momentum().norm())
      , q(par.charge())
      , navDir(ndir)
      , stepSize(ndir * std::abs(ssize))
    {
      // remember the start parameters
      startPos = pos;
      startDir = dir;
      // Init the jacobian matrix if needed
      if (par.covariance()) {
        // Get the reference surface for navigation
        const auto& surface = par.referenceSurface();
        // set the covariance transport flag to true and copy
        covTransport = true;
        cov          = ActsSymMatrixD<5>(*par.covariance());
        surface.initJacobianToGlobal(jacToGlobal, pos, dir, par.parameters());
      }
    }

    /// Global particle position accessor
    Vector3D
    position() const
    {
      return pos;
    }

    /// Momentum direction accessor
    Vector3D
    direction() const
    {
      return dir;
    }

    /// Actual momentum accessor
    Vector3D
    momentum() const
    {
      return p * dir;
    }

    /// Charge access
    double
    charge() const
    {
      return q;
    }

    /// Return a corrector
    corrector_t
    corrector() const
    {
      return corrector_t(startPos, startDir, pathAccumulated);
    }

    /// Method to update momentum, direction and p
    ///
    /// @param uposition the updated position
    /// @param udirection the updated direction
    /// @param up the updated momentum value
    void
    update(const Vector3D& uposition, const Vector3D& udirection, double up)
    {
      pos = uposition;
      dir = udirection;
      p   = up;
    }

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the state
    ///
    /// @param reinitialize is a flag to steer whether the
    ///        state should be reinitialized at the new
    ///        position
    ///
    /// @return the full transport jacobian
    const ActsMatrixD<5, 5>
    covarianceTransport(bool reinitialize = false)
    {
      // Optimized trigonometry on the propagation direction
      const double x = dir(0);  // == cos(phi) * sin(theta)
      const double y = dir(1);  // == sin(phi) * sin(theta)
      const double z = dir(2);  // == cos(theta)
      // can be turned into cosine/sine
      const double cosTheta    = z;
      const double sinTheta    = sqrt(x * x + y * y);
      const double invSinTheta = 1. / sinTheta;
      const double cosPhi      = x * invSinTheta;
      const double sinPhi      = y * invSinTheta;
      // prepare the jacobian to curvilinear
      ActsMatrixD<5, 7> jacToCurv = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
        // We normally operate in curvilinear coordinates defined as follows
        jacToCurv(0, 0) = -sinPhi;
        jacToCurv(0, 1) = cosPhi;
        jacToCurv(1, 0) = -cosPhi * cosTheta;
        jacToCurv(1, 1) = -sinPhi * cosTheta;
        jacToCurv(1, 2) = sinTheta;
      } else {
        // Under grazing incidence to z, the above coordinate system definition
        // becomes numerically unstable, and we need to switch to another one
        const double c    = sqrt(y * y + z * z);
        const double invC = 1. / c;
        jacToCurv(0, 1) = -z * invC;
        jacToCurv(0, 2) = y * invC;
        jacToCurv(1, 0) = c;
        jacToCurv(1, 1) = -x * y * invC;
        jacToCurv(1, 2) = -x * z * invC;
      }
      // Directional and momentum parameters for curvilinear
      jacToCurv(2, 3) = -sinPhi * invSinTheta;
      jacToCurv(2, 4) = cosPhi * invSinTheta;
      jacToCurv(3, 5) = -invSinTheta;
      jacToCurv(4, 6) = 1;
      // Apply the transport from the steps on the jacobian
      jacToGlobal = jacTransport * jacToGlobal;
      // Transport the covariance
      ActsRowVectorD<3>       normVec(dir);
      const ActsRowVectorD<5> sfactors
          = normVec * jacToGlobal.topLeftCorner<3, 5>();
      // The full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jacFull
          = jacToCurv * (jacToGlobal - derivative * sfactors);
      // Apply the actual covariance transport
      cov = (jacFull * cov * jacFull.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jacToGlobal  = ActsMatrixD<7, 5>::Zero();
        jacTransport = ActsMatrixD<7, 7>::Identity();
        // fill the jacobian to global for next transport
        jacToGlobal(0, eLOC_0) = -sinPhi;
        jacToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
        jacToGlobal(1, eLOC_0) = cosPhi;
        jacToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
        jacToGlobal(2, eLOC_1) = sinTheta;
        jacToGlobal(3, ePHI)   = -sinTheta * sinPhi;
        jacToGlobal(3, eTHETA) = cosTheta * cosPhi;
        jacToGlobal(4, ePHI)   = sinTheta * cosPhi;
        jacToGlobal(4, eTHETA) = cosTheta * sinPhi;
        jacToGlobal(5, eTHETA) = -sinTheta;
        jacToGlobal(6, eQOP)   = 1;
      }
      // return the full transport jacobian
      return jacFull;
    }

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the state
    ///
    /// @tparam S the Surfac type
    ///
    /// @param surface is the surface to which the covariance is
    ///        forwarded to
    /// @param reinitialize is a flag to steer whether the
    ///        state should be reinitialized at the new
    ///        position
    /// @note no check is done if the position is actually on the surface
    ///
    /// @return the full transport jacobian
    template <typename S>
    const ActsMatrixD<5, 5>
    covarianceTransport(const S& surface, bool reinitialize = false)
    {
      using VectorHelpers::phi;
      using VectorHelpers::theta;
      // Initialize the transport final frame jacobian
      ActsMatrixD<5, 7> jacToLocal = ActsMatrixD<5, 7>::Zero();
      // initalize the jacobian to local, returns the transposed ref frame
      auto rframeT = surface.initJacobianToLocal(jacToLocal, pos, dir);
      // Update the jacobian with the transport from the steps
      jacToGlobal = jacTransport * jacToGlobal;
      // calculate the form factors for the derivatives
      const ActsRowVectorD<5> sVec
          = surface.derivativeFactors(pos, dir, rframeT, jacToGlobal);
      // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jacFull
          = jacToLocal * (jacToGlobal - derivative * sVec);
      // Apply the actual covariance transport
      cov = (jacFull * cov * jacFull.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jacToGlobal  = ActsMatrixD<7, 5>::Zero();
        jacTransport = ActsMatrixD<7, 7>::Identity();
        // reset the derivative
        derivative = ActsVectorD<7>::Zero();
        // fill the jacobian to global for next transport
        Vector2D loc{0., 0.};
        surface.globalToLocal(pos, dir, loc);
        ActsVectorD<5> pars;
        pars << loc[eLOC_0], loc[eLOC_1], phi(dir), theta(dir), q / p;
        surface.initJacobianToGlobal(jacToGlobal, pos, dir, pars);
      }
      // store in the global jacobian
      jacobian = jacFull * jacobian;
      // return the full transport jacobian
      return jacFull;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0., 0., 0.);
    /// Global start particle position
    Vector3D startPos = Vector3D(0., 0., 0.);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1., 0., 0.);
    /// Momentum start direction (normalized)
    Vector3D startDir = Vector3D(1., 0., 0.);

    /// Momentum
    double p = 0.;

    /// The charge
    double q = 1.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// The full jacobian of the transport
    ActsMatrixD<5, 5> jacobian = ActsMatrixD<5, 5>::Identity();

    /// Jacobian from local to the global frame
    ActsMatrixD<7, 5> jacToGlobal = ActsMatrixD<7, 5>::Zero();

    /// Pure transport jacobian part from runge kutta integration
    ActsMatrixD<7, 7> jacTransport = ActsMatrixD<7, 7>::Identity();

    /// The propagation derivative
    ActsVectorD<7> derivative = ActsVectorD<7>::Zero();

    /// Covariance matrix (and indicator)
    //// associated with the initial error on track parameters
    bool              covTransport = false;
    ActsSymMatrixD<5> cov          = ActsSymMatrixD<5>::Zero();

    /// Lazily initialized state of the field Cache

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    typename BField::Cache fieldCache;

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
    cstep stepSize = std::numeric_limits<double>::max();

    /// List of algorithmic extensions
    extensionlist_t extension;

    /// Auctioneer for choosing the extension
    auctioneer_t auctioneer;
  };

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using state_type = State;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we usually return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor requires knowledge of the detector's magnetic field
  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  /// Convert the propagation state (global) to curvilinear parameters
  /// @param state The stepper state
  /// @param reinitialize is a flag to (optionally) reinitialse the state
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(State& state, bool reinitialize = false)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    // only do the covariance transport if needed
    if (state.covTransport) {
      // transport the covariance forward
      state.covarianceTransport(reinitialize);
      covPtr = std::make_unique<const ActsMatrixD<5, 5>>(state.cov);
    }
    // return the parameters
    return CurvilinearParameters(
        std::move(covPtr), state.pos, state.p * state.dir, state.q);
  }

  /// Convert the propagation state to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] state Propagation state used
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename S>
  static BoundParameters
  convert(State& state, const S& surface, bool reinitialize = false)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    // Perform error propagation if an initial covariance matrix was provided
    if (state.covTransport) {
      // transport the covariance forward
      state.covarianceTransport(surface, reinitialize);
      covPtr = std::make_unique<const ActsSymMatrixD<5>>(state.cov);
    }
    // return the bound parameters
    return BoundParameters(
        std::move(covPtr), state.pos, state.p * state.dir, state.q, surface);
  }

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(State& state, const Vector3D& pos) const
  {
    // get the field from the cell
    return m_bField.getField(pos, state.fieldCache);
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param[in,out] state is the propagation state associated with the track
  ///                      parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  double
  step(propagator_state_t& state) const
  {
    // Runge-Kutta integrator state
    StepData sd;

    double h2, half_h, error_estimate;

    // First Runge-Kutta point (at current position)
    sd.B_first = getField(state.stepping, state.stepping.pos);
    if (!state.stepping.extension.validExtensionForStep(state)
        || !state.stepping.extension.k1(state, sd.k1, sd.B_first)) {
      return 0.;
    }

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the local
    // integration error. The results are stated in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory
    const auto tryRungeKuttaStep = [&](const double h) -> bool {

      // State the square and half of the step size
      h2     = h * h;
      half_h = h * 0.5;

      // Second Runge-Kutta point
      const Vector3D pos1 = state.stepping.pos + half_h * state.stepping.dir
          + h2 * 0.125 * sd.k1;
      sd.B_middle = getField(state.stepping, pos1);
      if (!state.stepping.extension.k2(
              state, sd.k2, sd.B_middle, half_h, sd.k1)) {
        return false;
      }

      // Third Runge-Kutta point
      if (!state.stepping.extension.k3(
              state, sd.k3, sd.B_middle, half_h, sd.k2)) {
        return false;
      }

      // Last Runge-Kutta point
      const Vector3D pos2
          = state.stepping.pos + h * state.stepping.dir + h2 * 0.5 * sd.k3;
      sd.B_last = getField(state.stepping, pos2);
      if (!state.stepping.extension.k4(state, sd.k4, sd.B_last, h, sd.k3)) {
        return false;
      }

      // Return an estimate of the local integration error
      error_estimate = std::max(
          h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>(), 1e-20);
      return true;
    };

    // Select and adjust the appropriate Runge-Kutta step size as given
    // ATL-SOFT-PUB-2009-001
    while (!tryRungeKuttaStep(state.stepping.stepSize)
           || error_estimate > state.options.tolerance) {
      state.stepping.stepSize = state.stepping.stepSize
          * std::min(std::max(0.25,
                              std::pow((state.options.tolerance
                                        / std::abs(error_estimate)),
                                       0.25)),
                     4.);
      // If step size becomes too small the particle remains at the initial
      // place
      if (state.stepping.stepSize < state.options.stepSizeCutOff) {
        return 0.;  // Not moving due to too low momentum needs an aborter
      }
    }

    // use the adjusted step size
    const double h = state.stepping.stepSize;

    // When doing error propagation, update the associated Jacobian matrix
    if (state.stepping.covTransport) {
      // The step transport matrix in global coordinates
      ActsMatrixD<7, 7> D;
      if (!state.stepping.extension.finalize(state, h, sd, D)) {
        return 0.;
      }

      // for moment, only update the transport part
      state.stepping.jacTransport = D * state.stepping.jacTransport;
    } else {
      if (!state.stepping.extension.finalize(state, h)) {
        return 0.;
      }
    }

    // Update the track parameters according to the equations of motion
    state.stepping.pos
        += h * state.stepping.dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
    state.stepping.dir += h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
    state.stepping.dir /= state.stepping.dir.norm();
    state.stepping.derivative.template head<3>()     = state.stepping.dir;
    state.stepping.derivative.template segment<3>(3) = sd.k4;
    state.stepping.pathAccumulated += h;

    return h;
  }

private:
  /// Magnetic field inside of the detector
  BField m_bField;
};
}  // namespace Acts
