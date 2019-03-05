// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// LineSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

#ifdef ACTS_LINESURFACE_TRANSFORMATION_PLUGIN
#include ACTS_LINESURFACE_TRANSFORMATION_PLUGIN
#endif

inline void
LineSurface::localToGlobal(const GeometryContext& gctx,
                           const Vector2D&        lpos,
                           const Vector3D&        mom,
                           Vector3D&              gpos) const
{

  const auto& sTransform = transform(gctx);
  const auto& tMatrix    = sTransform.matrix();
  Vector3D    lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));

  // get the vector perpendicular to the momentum and the straw axis
  Vector3D radiusAxisGlobal(lineDirection.cross(mom));
  Vector3D locZinGlobal = sTransform * Vector3D(0., 0., lpos[eLOC_Z]);
  // add eLOC_R * radiusAxis
  gpos = Vector3D(locZinGlobal + lpos[eLOC_R] * radiusAxisGlobal.normalized());
}

inline bool
LineSurface::globalToLocal(const GeometryContext& gctx,
                           const Vector3D&        gpos,
                           const Vector3D&        mom,
                           Vector2D&              lpos) const
{
  using VectorHelpers::perp;

  const auto& sTransform = transform(gctx);
  const auto& tMatrix    = sTransform.matrix();
  Vector3D    lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  // Bring the global position into the local frame
  Vector3D loc3Dframe = sTransform.inverse() * gpos;
  // construct localPosition with sign*perp(candidate) and z.()
  lpos = Vector2D(perp(loc3Dframe), loc3Dframe.z());
  Vector3D sCenter(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
  Vector3D decVec(gpos - sCenter);
  // assign the right sign
  double sign = ((lineDirection.cross(mom)).dot(decVec) < 0.) ? -1. : 1.;
  lpos[eLOC_R] *= sign;
  return true;
}

inline std::string
LineSurface::name() const
{
  return "Acts::LineSurface";
}

inline const RotationMatrix3D
LineSurface::referenceFrame(const GeometryContext& gctx,
                            const Vector3D& /*unused*/,
                            const Vector3D& mom) const
{
  RotationMatrix3D mFrame;
  const auto&      tMatrix = transform(gctx).matrix();
  Vector3D         measY(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  Vector3D         measX(measY.cross(mom).normalized());
  Vector3D         measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

inline double
LineSurface::pathCorrection(const GeometryContext& /*unused*/,
                            const Vector3D& /*pos*/,
                            const Vector3D& /*mom*/) const
{
  return 1.;
}

inline const Vector3D
LineSurface::binningPosition(const GeometryContext& gctx,
                             BinningValue /*bValue*/) const
{
  return center(gctx);
}

inline const Vector3D
LineSurface::normal(const GeometryContext& gctx, const Vector2D& /*lpos*/) const
{
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const SurfaceBounds&
LineSurface::bounds() const
{
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

inline Intersection
LineSurface::intersectionEstimate(const GeometryContext& gctx,
                                  const Vector3D&        gpos,
                                  const Vector3D&        gdir,
                                  NavigationDirection    navDir,
                                  const BoundaryCheck&   bcheck,
                                  CorrFnc                correct) const
{
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  Vector3D ma = gpos;
  Vector3D ea = gdir;
  // line two is the line surface
  const auto& tMatrix = transform(gctx).matrix();
  Vector3D    mb      = tMatrix.block<3, 1>(0, 3).transpose();
  Vector3D    eb      = tMatrix.block<3, 1>(0, 2).transpose();
  // now go ahead and solve for the closest approach
  Vector3D mab(mb - ma);
  double   eaTeb = ea.dot(eb);
  double   denom = 1 - eaTeb * eaTeb;
  // validity parameter
  bool valid = false;
  if (denom * denom > s_onSurfaceTolerance * s_onSurfaceTolerance) {
    double u = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // evaluate in terms of direction
    valid = (navDir * u >= 0);
    // evaluate validaty in terms of bounds
    Vector3D result = (ma + u * ea);
    // update if you have a correction
    if (correct && correct(ma, ea, u)) {
      // update everything that is in relation to ea
      eaTeb = ea.dot(eb);
      denom = 1 - eaTeb * eaTeb;
      if (denom * denom > s_onSurfaceTolerance * s_onSurfaceTolerance) {
        u      = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
        result = (ma + u * ea);
        // if you have specified a navigation direction, valid mean path > 0.
        valid = (navDir * u >= 0);
      } else {
        valid = false;
      }
    }
    // it just needs to be a insideBounds() check
    // @todo there should be a faster check possible
    valid = bcheck ? (valid && isOnSurface(gctx, result, gdir, bcheck)) : valid;
    // return the result with validity
    return Intersection(result, u, valid);
  }
  // return a false intersection
  return Intersection(gpos, std::numeric_limits<double>::max(), false);
}

//~ inline void LineSurface::initJacobianToGlobal(const GeometryContext& gctx, TrackToGlobalMatrix& jacobian,
                                              //~ const Vector3D&       gpos,
                                              //~ const Vector3D&       dir,
                                              //~ const TrackVector& pars) const

//~ inline const TrackRowVector
//~ LineSurface::derivativeFactors(const GeometryContext&  gctx,const Vector3D&         pos,
                               //~ const Vector3D&         dir,
                               //~ const RotationMatrix3D& rft,
                               //~ const TrackToGlobalMatrix& jac) const