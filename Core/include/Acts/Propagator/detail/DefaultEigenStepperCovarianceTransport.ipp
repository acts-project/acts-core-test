// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

void
covarianceTransport(State& state, bool reinitialize = false) const
{
  // Optimized trigonometry on the propagation direction
  const double x = state.dir(0);  // == cos(phi) * sin(theta)
  const double y = state.dir(1);  // == sin(phi) * sin(theta)
  const double z = state.dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta    = z;
  const double sinTheta    = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi      = x * invSinTheta;
  const double sinPhi      = y * invSinTheta;
  // prepare the jacobian to curvilinear
  GlobalToTrackMatrix jacToCurv = GlobalToTrackMatrix::Zero();
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
  state.jacToGlobal = state.jacTransport * state.jacToGlobal;
  // Transport the covariance
  ActsRowVectorD<3>    normVec(state.dir);
  const TrackRowVector sfactors
      = normVec * state.jacToGlobal.template topLeftCorner<3, TrackParsDim>();
  // The full jacobian is ([to local] jacobian) * ([transport] jacobian)
  const TrackMatrix jacFull
      = jacToCurv * (state.jacToGlobal - state.derivative * sfactors);
  // Apply the actual covariance transport
  state.cov = (jacFull * state.cov * jacFull.transpose());
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    // reset the jacobians
    state.jacToGlobal  = TrackToGlobalMatrix::Zero();
    state.jacTransport = GlobalMatrix::Identity();
    // fill the jacobian to global for next transport
    state.jacToGlobal(0, eLOC_0) = -sinPhi;
    state.jacToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
    state.jacToGlobal(1, eLOC_0) = cosPhi;
    state.jacToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
    state.jacToGlobal(2, eLOC_1) = sinTheta;
    state.jacToGlobal(3, ePHI)   = -sinTheta * sinPhi;
    state.jacToGlobal(3, eTHETA) = cosTheta * cosPhi;
    state.jacToGlobal(4, ePHI)   = sinTheta * cosPhi;
    state.jacToGlobal(4, eTHETA) = cosTheta * sinPhi;
    state.jacToGlobal(5, eTHETA) = -sinTheta;
    state.jacToGlobal(6, eQOP)   = 1;
  }
  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull * state.jacobian;
}

template <typename surface_t>
void
covarianceTransport(State&           state,
                    const surface_t& surface,
                    bool             reinitialize = true) const
{
  using VectorHelpers::phi;
  using VectorHelpers::theta;
  // Initialize the transport final frame jacobian
  GlobalToTrackMatrix jacToLocal = GlobalToTrackMatrix::Zero();
  // initalize the jacobian to local, returns the transposed ref frame
  auto rframeT = surface.initJacobianToLocal(jacToLocal, state.pos, state.dir);
  // Update the jacobian with the transport from the steps
  state.jacToGlobal = state.jacTransport * state.jacToGlobal;
  // calculate the form factors for the derivatives
  const TrackRowVector sVec = surface.derivativeFactors(
      state.pos, state.dir, rframeT, state.jacToGlobal);
  // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
  const TrackMatrix jacFull
      = jacToLocal * (state.jacToGlobal - state.derivative * sVec);
  // Apply the actual covariance transport
  state.cov = (jacFull * state.cov * jacFull.transpose());
  state.cov(5, 5) = 1;  // TODO: get rid of this; it just allows the
                        // calculation of inv(cov) for the GainMatrixSmoother
                        // variable G
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    // reset the jacobians
    state.jacToGlobal  = TrackToGlobalMatrix::Zero();
    state.jacTransport = GlobalMatrix::Identity();
    // reset the derivative
    state.derivative = GlobalVector::Zero();
    // fill the jacobian to global for next transport
    Vector2D loc{0., 0.};
    surface.globalToLocal(state.pos, state.dir, loc);
    TrackVector pars = TrackVector::Zero();
    pars(0)          = loc[eLOC_0];
    pars(1)          = loc[eLOC_1];
    pars(2)          = phi(state.dir);
    pars(3)          = theta(state.dir);
    pars(4)          = state.q / state.p;
    surface.initJacobianToGlobal(state.jacToGlobal, state.pos, state.dir, pars);
  }
  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull * state.jacobian;
}