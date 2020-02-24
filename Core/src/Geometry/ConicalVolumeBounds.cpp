// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ConicalVolumeBounds.hpp"
#include <cmath>
#include <iostream>
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

Acts::ConicalVolumeBounds::ConicalVolumeBounds(
    double innerAlpha, double innerOffsetZ, double outerAlpha,
    double outerOffsetZ, double halflengthZ, double averagePhi,
    double halfPhiSector)
    : VolumeBounds(), m_valueStore(bv_length, 0.) {
  VOLUMEBOUNDS_VALUESTORE_FILL(innerAlpha);
  VOLUMEBOUNDS_VALUESTORE_FILL(innerOffsetZ);
  VOLUMEBOUNDS_VALUESTORE_FILL(outerAlpha);
  VOLUMEBOUNDS_VALUESTORE_FILL(outerOffsetZ);
  VOLUMEBOUNDS_VALUESTORE_FILL(halflengthZ);
  VOLUMEBOUNDS_VALUESTORE_FILL(averagePhi);
  VOLUMEBOUNDS_VALUESTORE_FILL(halfPhiSector);

  if (innerAlpha > 0.) {
    m_innerTanAlpha = std::tan(innerAlpha);
    double innerZmin = innerOffsetZ - halflengthZ;
    double innerZmax = innerOffsetZ + halflengthZ;
    m_innerRmin = std::abs(innerZmin) * m_innerTanAlpha;
    m_innerRmax = std::abs(innerZmax) * m_innerTanAlpha;
    m_innerConeBounds = std::make_shared<ConeBounds>(
        innerAlpha, innerZmin, innerZmax, halfPhiSector, averagePhi);
  }

  m_outerTanAlpha = std::tan(outerAlpha);
  double outerZmin = outerOffsetZ - halflengthZ;
  double outerZmax = outerOffsetZ + halflengthZ;
  m_outerRmin = std::abs(outerZmin) * m_outerTanAlpha;
  m_outerRmax = std::abs(outerZmax) * m_outerTanAlpha;

  if (halflengthZ < outerOffsetZ) {
    m_negativeDiscBounds = std::make_shared<RadialBounds>(
        innerRmin(), outerRmin(), averagePhi, halfPhiSector);
  }

  m_positiveDiscBounds = std::make_shared<RadialBounds>(
      innerRmax(), outerRmax(), averagePhi, halfPhiSector);
  m_outerConeBounds = std::make_shared<ConeBounds>(
      outerAlpha, outerZmin, outerZmax, halfPhiSector, averagePhi);

  createSectorBounds();
}

Acts::ConicalVolumeBounds::ConicalVolumeBounds(double cylinderR, double alpha,
                                               double offsetZ,
                                               double halflengthZ,
                                               double averagePhi,
                                               double halfPhiSector)
    : VolumeBounds(), m_valueStore(bv_length, 0.) {
  VOLUMEBOUNDS_VALUESTORE_FILL(halflengthZ);
  VOLUMEBOUNDS_VALUESTORE_FILL(averagePhi);
  VOLUMEBOUNDS_VALUESTORE_FILL(halfPhiSector);

  // Cone parameters
  double tanAlpha = std::tan(alpha);
  double zmin = offsetZ - halflengthZ;
  double zmax = offsetZ + halflengthZ;
  double rmin = std::abs(zmin) * tanAlpha;
  double rmax = std::abs(zmax) * tanAlpha;

  if (rmin >= cylinderR) {
    // Cylindrical cut-out of a cone
    m_innerRmin = cylinderR;
    m_innerRmax = cylinderR;
    m_innerCylinderBounds = std::make_shared<CylinderBounds>(
        cylinderR, averagePhi, halfPhiSector, halflengthZ);
    m_outerTanAlpha = tanAlpha;
    m_outerRmin = rmin;
    m_outerRmax = rmax;
    m_valueStore[bv_outerOffsetZ] = offsetZ;
    m_outerConeBounds = std::make_shared<ConeBounds>(alpha, zmin, zmax,
                                                     halfPhiSector, averagePhi);
  } else if (rmax <= cylinderR) {
    // Conical cut-out of a cylinder
    m_outerRmin = cylinderR;
    m_outerRmax = cylinderR;
    m_outerCylinderBounds = std::make_shared<CylinderBounds>(
        cylinderR, averagePhi, halfPhiSector, halflengthZ);
    m_innerTanAlpha = tanAlpha;
    m_innerRmin = rmin;
    m_innerRmax = rmax;
    m_valueStore[bv_innerOffsetZ] = offsetZ;
    m_innerConeBounds = std::make_shared<ConeBounds>(alpha, zmin, zmax,
                                                     halfPhiSector, averagePhi);
  } else {
    throw std::domain_error(
        "Cylinder and Cone are intersecting, not possible.");
  }

  if (halflengthZ < offsetZ) {
    m_negativeDiscBounds = std::make_shared<RadialBounds>(
        innerRmin(), outerRmin(), averagePhi, halfPhiSector);
  }

  m_positiveDiscBounds = std::make_shared<RadialBounds>(
      innerRmax(), outerRmax(), averagePhi, halfPhiSector);

  createSectorBounds();
}

Acts::ConicalVolumeBounds::ConicalVolumeBounds(const ConicalVolumeBounds& bobo)
    : VolumeBounds(), m_valueStore(bobo.m_valueStore) {}

Acts::ConicalVolumeBounds* Acts::ConicalVolumeBounds::clone() const {
  return new ConicalVolumeBounds(*this);
}

Acts::ConicalVolumeBounds& Acts::ConicalVolumeBounds::operator=(
    const ConicalVolumeBounds& cobo) {
  if (this != &cobo) {
    m_valueStore = cobo.m_valueStore;
  }
  return *this;
}

Acts::SurfacePtrVector Acts::ConicalVolumeBounds::decomposeToSurfaces(
    const Transform3D* transformPtr) const {
  // The transform - apply when given
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);

  SurfacePtrVector rSurfaces;
  rSurfaces.reserve(6);

  // Create an inner Cone
  if (m_innerConeBounds != nullptr) {
    auto innerConeRelTrans = Transform3D::Identity();
    innerConeRelTrans.pretranslate(Vector3D(0., 0., -innerOffsetZ()));
    auto innerConeAbsTrans =
        std::make_shared<Transform3D>(transform * innerConeRelTrans);
    auto innerCone =
        Surface::makeShared<ConeSurface>(innerConeAbsTrans, m_innerConeBounds);
    rSurfaces.push_back(innerCone);
  } else if (m_innerCylinderBounds != nullptr) {
    // or alternatively the inner Cylinder
    auto innerCylinderAbsTrans = std::make_shared<Transform3D>(transform);
    auto innerCylinder = Surface::makeShared<CylinderSurface>(
        innerCylinderAbsTrans, m_innerCylinderBounds);
    rSurfaces.push_back(innerCylinder);
  }

  // Create an outer Cone
  if (m_outerConeBounds != nullptr) {
    auto outerConeRelTrans = Transform3D::Identity();
    outerConeRelTrans.pretranslate(Vector3D(0., 0., -outerOffsetZ()));
    auto outerConeAbsTrans =
        std::make_shared<Transform3D>(transform * outerConeRelTrans);
    auto outerCone =
        Surface::makeShared<ConeSurface>(outerConeAbsTrans, m_outerConeBounds);
    rSurfaces.push_back(outerCone);
  } else if (m_outerCylinderBounds != nullptr) {
    // or alternatively an outer Cylinder
    auto outerCylinderAbsTrans = std::make_shared<Transform3D>(transform);
    auto outerCylinder = Surface::makeShared<CylinderSurface>(
        outerCylinderAbsTrans, m_outerCylinderBounds);
    rSurfaces.push_back(outerCylinder);
  }

  // Set a disc at Zmin
  if (m_negativeDiscBounds != nullptr) {
    auto negativeDiscRelTrans = Transform3D::Identity();
    negativeDiscRelTrans.pretranslate(Vector3D(0., 0., -halflengthZ()));
    auto negativeDiscAbsTrans =
        std::make_shared<Transform3D>(transform * negativeDiscRelTrans);
    auto negativeDisc = Surface::makeShared<DiscSurface>(negativeDiscAbsTrans,
                                                         m_negativeDiscBounds);
    rSurfaces.push_back(negativeDisc);
  }

  // Set a disc at Zmax
  auto positiveDiscRelTrans = Transform3D::Identity();
  positiveDiscRelTrans.pretranslate(Vector3D(0., 0., halflengthZ()));
  auto positiveDiscAbsTrans =
      std::make_shared<Transform3D>(transform * positiveDiscRelTrans);
  auto positiveDisc = Surface::makeShared<DiscSurface>(positiveDiscAbsTrans,
                                                       m_positiveDiscBounds);
  rSurfaces.push_back(positiveDisc);

  if (m_sectorBounds) {
    RotationMatrix3D sectorRotation;
    sectorRotation.col(0) = s_zAxis;
    sectorRotation.col(1) = s_xAxis;
    sectorRotation.col(2) = s_yAxis;

    // curvilinear surfaces are boundless
    Transform3D negSectorRelTrans{sectorRotation};
    negSectorRelTrans.prerotate(
        AngleAxis3D(averagePhi() - halfPhiSector(), s_zAxis));
    auto negSectorAbsTrans =
        std::make_shared<Transform3D>(transform * negSectorRelTrans);

    rSurfaces.push_back(
        Surface::makeShared<PlaneSurface>(negSectorAbsTrans, m_sectorBounds));

    Transform3D posSectorRelTrans{sectorRotation};
    posSectorRelTrans.prerotate(
        AngleAxis3D(averagePhi() + halfPhiSector(), s_zAxis));
    auto posSectorAbsTrans =
        std::make_shared<Transform3D>(transform * posSectorRelTrans);

    rSurfaces.push_back(
        Surface::makeShared<PlaneSurface>(posSectorAbsTrans, m_sectorBounds));
  }
  return rSurfaces;
}
bool Acts::ConicalVolumeBounds::inside(const Vector3D& pos, double tol) const {
  double z = pos.z();
  double zmin = z + tol;
  double zmax = z - tol;
  // Quick check ouside z
  if (zmin < -halflengthZ() or zmax > halflengthZ()) {
    return false;
  }
  double r = VectorHelpers::perp(pos);
  if (std::abs(halfPhiSector() - M_PI) > s_onSurfaceTolerance) {
    // need to check the phi sector - approximate phi tolerance
    double phitol = tol / r;
    double phi = VectorHelpers::phi(pos);
    double phimin = phi - phitol;
    double phimax = phi + phitol;
    if (phimin < averagePhi() - halfPhiSector() or
        phimax > averagePhi() + halfPhiSector()) {
      return false;
    }
  }
  // We are within phi sector check box r quickly
  double rmin = r + tol;
  double rmax = r - tol;
  if (rmin > innerRmax() and rmax < outerRmin()) {
    return true;
  }
  // Finally we need to check the cone
  if (m_innerConeBounds != nullptr) {
    double innerConeR = m_innerConeBounds->r(std::abs(z + innerOffsetZ()));
    if (innerConeR > rmin) {
      return false;
    }
  } else if (innerRmax() > rmin) {
    return false;
  }
  // And the outer cone
  if (m_outerConeBounds != nullptr) {
    double outerConeR = m_outerConeBounds->r(std::abs(z + outerOffsetZ()));
    if (outerConeR < rmax) {
      return false;
    }
  } else if (outerRmax() < rmax) {
    return false;
  }
  return true;
}

void Acts::ConicalVolumeBounds::createSectorBounds() {
  if (std::abs(halfPhiSector() - M_PI) > s_onSurfaceTolerance) {
    // The 4 points building the sector
    std::vector<Vector2D> polyVertices = {{-halflengthZ(), innerRmin()},
                                          {halflengthZ(), innerRmax()},
                                          {halflengthZ(), outerRmax()},
                                          {-halflengthZ(), outerRmin()}};
    m_sectorBounds =
        std::make_shared<ConvexPolygonBounds<4>>(std::move(polyVertices));
  }
}

// ostream operator overload
std::ostream& Acts::ConicalVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT(sl);
}

Acts::Volume::BoundingBox Acts::ConicalVolumeBounds::boundingBox(
    const Acts::Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  Vector3D vmin(-outerRmax(), -outerRmax(), -0.5 * halflengthZ());
  Vector3D vmax(outerRmax(), outerRmax(), 0.5 * halflengthZ());
  Volume::BoundingBox box(entity, vmin - envelope, vmax + envelope);
  return trf == nullptr ? box : box.transformed(*trf);
}
