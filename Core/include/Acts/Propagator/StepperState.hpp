// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>
#include <variant>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// State for track parameter propagation
///
struct StepperState {
  using Jacobian = std::variant<BoundMatrix, FreeToBoundMatrix,
                                BoundToFreeMatrix, FreeMatrix>;
  using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;

  /// Delete the default constructor
  StepperState() = delete;

  /// Constructor from the initial track parameters
  ///
  /// @tparam parameters_t the Type of the track parameters
  ///
  /// @param [in] gctx is the context object for the geometery
  /// @param [in] mctx is the context object for the magnetic field
  /// @param [in] par The track parameters at start
  /// @param [in] ndir is the navigation direction
  /// @param [in] ssize is the (absolute) maximum step size
  /// @param [in] stolerance is the stepping tolerance
  template <typename parameters_t,
            std::enable_if_t<parameters_t::is_local_representation, int> = 0>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max(),
      double stolerance = s_onSurfaceTolerance)
      : pos(par.position()),
        dir(par.momentum().normalized()),
        p(par.momentum().norm()),
        q(par.charge()),
        t(par.time()),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    if (par.covariance()) {
      // Set the covariance transport flag to true
      covTransport = true;
      // Get the covariance
      jacToGlobal = BoundToFreeMatrix::Zero();
      par.referenceSurface().initJacobianToGlobal(gctx, *jacToGlobal, pos, dir,
                                                  par.parameters());
      cov = *par.covariance();
      jacobian.emplace<0>(BoundMatrix::Identity());
    }
  }

  /// Constructor from the initial track parameters
  ///
  /// @tparam parameters_t the Type of the track parameters
  ///
  /// @param [in] gctx is the context object for the geometery
  /// @param [in] mctx is the context object for the magnetic field
  /// @param [in] par The track parameters at start
  /// @param [in] ndir is the navigation direction
  /// @param [in] ssize is the (absolute) maximum step size
  /// @param [in] stolerance is the stepping tolerance
  template <
      typename parameters_t,
      std::enable_if_t<not parameters_t::is_local_representation, int> = 0>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max(),
      double stolerance = s_onSurfaceTolerance)
      : pos(par.position()),
        dir(par.momentum().normalized()),
        p(par.momentum().norm()),
        q(par.charge()),
        t(par.time()),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    if (par.covariance()) {
      // Set the covariance transport flag to true
      covTransport = true;
      // Get the covariance
      cov = *par.covariance();
      jacobian.emplace<3>(FreeMatrix::Identity());
    }
  }

  /// Jacobian from local to the global frame
  std::optional<BoundToFreeMatrix> jacToGlobal;

  /// Pure transport jacobian part from runge kutta integration
  FreeMatrix jacTransport = FreeMatrix::Identity();

  /// The full jacobian since the first step
  Jacobian jacobian;

  /// The propagation derivative
  FreeVector derivative = FreeVector::Zero();

  /// Boolean to indiciate if you need covariance transport
  bool covTransport = false;
  Covariance cov;

  /// Global particle position
  Vector3D pos = Vector3D(0., 0., 0.);

  /// Momentum direction (normalized)
  Vector3D dir = Vector3D(1., 0., 0.);

  /// Momentum
  double p = 0.;

  /// Save the charge: neutral as default for SL stepper
  double q = 0.;

  /// Propagated time
  double t = 0.;

  /// Navigation direction, this is needed for searching
  NavigationDirection navDir;

  /// accummulated path length state
  double pathAccumulated = 0.;

  /// adaptive step size of the runge-kutta integration
  ConstrainedStep stepSize = std::numeric_limits<double>::max();
  /// Previous step size for overstep estimation (ignored for SL stepper)
  double previousStepSize = 0.;

  /// The tolerance for the stepping
  double tolerance = s_onSurfaceTolerance;

  /// Cache the geometry context of this propagation
  std::reference_wrapper<const GeometryContext> geoContext;
};
}  // namespace Acts