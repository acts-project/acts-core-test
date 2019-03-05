// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Surface.ipp, Acts project
///////////////////////////////////////////////////////////////////

#ifdef ACTS_SURFACE_TRANSFORMATION_PLUGIN
#include ACTS_SURFACE_TRANSFORMATION_PLUGIN
#endif

inline const Vector3D
Surface::center(const GeometryContext& gctx) const
{
  // fast access via tranform matrix (and not translation())
  auto tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
}

inline const Acts::Vector3D
Surface::normal(const GeometryContext& gctx, const Vector3D& /*unused*/) const
{
  return normal(gctx, s_origin2D);
}

inline const Transform3D&
Surface::transform(const GeometryContext& gctx) const
{
  if (m_transform != nullptr) {
    return (*(m_transform.get()));
  }
  if (m_associatedDetElement != nullptr) {
    return m_associatedDetElement->transform(gctx);
  }
  return s_idTransform;
}

inline bool
Surface::insideBounds(const Vector2D& locpos, const BoundaryCheck& bcheck) const
{
  return bounds().inside(locpos, bcheck);
}

inline const RotationMatrix3D
Surface::referenceFrame(const GeometryContext& gctx,
                        const Vector3D& /*unused*/,
                        const Vector3D& /*unused*/) const
{
  return transform(gctx).matrix().block<3, 3>(0, 0);
}

//~ inline void Surface::initJacobianToGlobal(const GeometryContext& gctx, TrackToGlobalMatrix& jacobian,
                                          //~ const Vector3D& gpos,
                                          //~ const Vector3D& dir,
                                          //~ const TrackVector& /*pars*/) const

//~ inline const RotationMatrix3D
    //~ Surface::initJacobianToLocal(const GeometryContext& gctx, GlobalToTrackMatrix& jacobian,
                                 //~ const Vector3D& gpos,
                                 //~ const Vector3D& dir) const

//~ inline const TrackRowVector
//~ Surface::derivativeFactors(const GeometryContext& /*unused*/, const Vector3D& /*unused*/,
                           //~ const Vector3D&         dir,
                           //~ const RotationMatrix3D& rft,
                           //~ const TrackToGlobalMatrix& jac) const

template <typename parameters_t>
bool
Surface::isOnSurface(const GeometryContext& gctx,
                     const parameters_t&    pars,
                     const BoundaryCheck&   bcheck) const
{
  // surface pointer comparison as a first fast check (w/o transform)
  // @todo check if we can find a fast way that works for stepper state and
  // parameters
  // if ((&pars.referenceSurface() == this) && !bcheck) return true;
  return isOnSurface(gctx, pars.position(), pars.momentum(), bcheck);
}

inline const DetectorElementBase*
Surface::associatedDetectorElement() const
{
  return m_associatedDetElement;
}

inline const Layer*
Surface::associatedLayer() const
{
  return (m_associatedLayer);
}

inline const ISurfaceMaterial*
Surface::surfaceMaterial() const
{
  return m_surfaceMaterial.get();
}

inline void
Surface::assignSurfaceMaterial(std::shared_ptr<const ISurfaceMaterial> material)
{
  m_surfaceMaterial = std::move(material);
}

inline void
Surface::associateLayer(const Layer& lay)
{
  m_associatedLayer = (&lay);
}
