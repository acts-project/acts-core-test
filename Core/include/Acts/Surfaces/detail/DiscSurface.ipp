// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// DiscSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

#ifdef ACTS_DISCSURFACE_TRANSFORMATION_PLUGIN
#include ACTS_DISCSURFACE_TRANSFORMATION_PLUGIN
#endif

inline const Vector2D
DiscSurface::localPolarToCartesian(const Vector2D& lpolar) const
{
  return Vector2D(lpolar[eLOC_R] * cos(lpolar[eLOC_PHI]),
                  lpolar[eLOC_R] * sin(lpolar[eLOC_PHI]));
}

inline const Vector2D
DiscSurface::localCartesianToPolar(const Vector2D& lcart) const
{
  return Vector2D(
      sqrt(lcart[eLOC_X] * lcart[eLOC_X] + lcart[eLOC_Y] * lcart[eLOC_Y]),
      atan2(lcart[eLOC_Y], lcart[eLOC_X]));
}

//~ inline void DiscSurface::initJacobianToGlobal(const GeometryContext& gctx, TrackToGlobalMatrix& jacobian,
                                              //~ const Vector3D&       gpos,
                                              //~ const Vector3D&       dir,
                                              //~ const TrackVector& pars) const


//~ inline const RotationMatrix3D
    //~ DiscSurface::initJacobianToLocal(const GeometryContext& gctx, GlobalToTrackMatrix& jacobian,
                                     //~ const Vector3D& gpos,
                                     //~ const Vector3D& dir) const

inline Intersection
DiscSurface::intersectionEstimate(const GeometryContext& gctx,
                                  const Vector3D&        gpos,
                                  const Vector3D&        gdir,
                                  NavigationDirection    navDir,
                                  const BoundaryCheck&   bcheck,
                                  CorrFnc                correct) const
{
  // minimize the call to transform()
  const auto&    tMatrix = transform(gctx).matrix();
  const Vector3D pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3D pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // return solution and path
  Vector3D solution(0., 0., 0.);
  double   path = std::numeric_limits<double>::infinity();
  // lemma : the solver -> should catch current values
  auto solve = [&solution, &path, &pnormal, &pcenter, &navDir](
      const Vector3D& lpos, const Vector3D& ldir) -> bool {
    double denom = ldir.dot(pnormal);
    if (denom != 0.0) {
      path     = (pnormal.dot((pcenter - lpos))) / (denom);
      solution = (lpos + path * ldir);
    }
    // is valid if it goes into the right direction
    return ((navDir == 0) || path * navDir >= 0.);
  };
  // solve first
  bool valid = solve(gpos, gdir);
  // if configured to correct, do it and solve again
  if (correct) {
    // copy as the corrector may change them
    Vector3D lposc = gpos;
    Vector3D ldirc = gdir;
    if (correct(lposc, ldirc, path)) {
      valid = solve(lposc, ldirc);
    }
  }
  // evaluate (if necessary in terms of boundaries)
  // @todo: speed up isOnSurface - we know that it is on surface
  //  all we need is to check if it's inside bounds in 3D space
  valid = bcheck ? (valid && isOnSurface(gctx, solution, gdir, bcheck)) : valid;
  // return the result
  return Intersection(solution, path, valid);
}

inline const Vector3D
DiscSurface::normal(const GeometryContext& gctx,
                    const Vector2D& /*unused*/) const
{
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D
DiscSurface::binningPosition(const GeometryContext& gctx,
                             BinningValue /*unused*/) const
{
  return center(gctx);
}

inline double
DiscSurface::pathCorrection(const GeometryContext& gctx,
                            const Vector3D&        pos,
                            const Vector3D&        mom) const
{
  /// we can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, pos).dot(mom.normalized()));
}
