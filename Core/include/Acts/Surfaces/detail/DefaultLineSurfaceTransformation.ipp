// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline void
LineSurface::initJacobianToGlobal(TrackToGlobalMatrix& jacobian,
                                  const Vector3D&      gpos,
                                  const Vector3D&      dir,
                                  const TrackVector&   pars) const
{
  // The trigonometry required to convert the direction to spherical
  // coordinates and then compute the sines and cosines again can be
  // surprisingly expensive from a performance point of view.
  //
  // Here, we can avoid it because the direction is by definition a unit
  // vector, with the following coordinate conversions...
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  const double z = dir(2);  // == cos(theta)

  // ...which we can invert to directly get the sines and cosines:
  const double cos_theta     = z;
  const double sin_theta     = sqrt(x * x + y * y);
  const double inv_sin_theta = 1. / sin_theta;
  const double cos_phi       = x * inv_sin_theta;
  const double sin_phi       = y * inv_sin_theta;
  // retrieve the reference frame
  const auto rframe = referenceFrame(gpos, dir);
  // the local error components - given by the reference frame
  jacobian.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the momentum components
  jacobian(3, ePHI)   = (-sin_theta) * sin_phi;
  jacobian(3, eTHETA) = cos_theta * cos_phi;
  jacobian(4, ePHI)   = sin_theta * cos_phi;
  jacobian(4, eTHETA) = cos_theta * sin_phi;
  jacobian(5, eTHETA) = (-sin_theta);
  jacobian(6, eQOP)   = 1;

  // the projection of direction onto ref frame normal
  double ipdn = 1. / dir.dot(rframe.col(2));
  // build the cross product of d(D)/d(ePHI) components with y axis
  auto dDPhiY = rframe.block<3, 1>(0, 1).cross(jacobian.block<3, 1>(3, ePHI));
  // and the same for the d(D)/d(eTheta) components
  auto dDThetaY
      = rframe.block<3, 1>(0, 1).cross(jacobian.block<3, 1>(3, eTHETA));
  // and correct for the x axis components
  dDPhiY -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDPhiY));
  dDThetaY
      -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDThetaY));
  // set the jacobian components for global d/ phi/Theta
  jacobian.block<3, 1>(0, ePHI)   = dDPhiY * pars[eLOC_0] * ipdn;
  jacobian.block<3, 1>(0, eTHETA) = dDThetaY * pars[eLOC_0] * ipdn;
}

inline const TrackRowVector
LineSurface::derivativeFactors(const Vector3D&            pos,
                               const Vector3D&            dir,
                               const RotationMatrix3D&    rft,
                               const TrackToGlobalMatrix& jac) const
{
  // the vector between position and center
  ActsRowVectorD<3> pc = (pos - center()).transpose();
  // the longitudinal component vector (alogn local z)
  ActsRowVectorD<3> locz = rft.block<1, 3>(1, 0);
  // build the norm vector comonent by subtracting the longitudinal one
  double            long_c   = locz * dir;
  ActsRowVectorD<3> norm_vec = dir.transpose() - long_c * locz;
  // calculate the s factors for the dependency on X
  const TrackRowVector s_vec = norm_vec * jac.topLeftCorner<3, TrackParsDim>();
  // calculate the d factors for the dependency on Tx
  const TrackRowVector d_vec = locz * jac.block<3, TrackParsDim>(3, 0);
  // normalisation of normal & longitudinal components
  double norm = 1. / (1. - long_c * long_c);
  // create a matrix representation
  ActsMatrixD<3, TrackParsDim> long_mat = ActsMatrixD<3, TrackParsDim>::Zero();
  long_mat.colwise() += locz.transpose();
  // build the combined normal & longitudinal components
  return (norm * (s_vec
                  - pc * (long_mat * d_vec.asDiagonal()
                          - jac.block<3, TrackParsDim>(3, 0))));
}
