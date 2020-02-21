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

  if (innerAlpha) {
    m_innerTanAlpha = std::tan(innerAlpha);
    double innerZmin = innerOffsetZ - halflengthZ;
    double innerZmax = innerOffsetZ + halflengthZ;
    m_innerRmin = std::abs(innerZmin) * m_innerTanAlpha;
    m_innerRmax = std::abs(innerZmax) * m_innerTanAlpha;
    m_innerConeBounds = std::make_shared<ConeBounds>(
        innerAlpha, averagePhi, innerZmin, innerZmax, halfPhiSector);
  }

  if (halflengthZ < outerOffsetZ) {
    m_smallerDiscBounds = std::make_shared<RadialBounds>(
        innerRmin(), outerRmin(), averagePhi, halfPhiSector);
  }

  m_outerTanAlpha = std::tan(outerAlpha);
  double outerZmin = outerOffsetZ - halflengthZ;
  double outerZmax = outerOffsetZ + halflengthZ;
  m_outerRmin = std::abs(outerZmin) * m_outerTanAlpha;
  m_outerRmax = std::abs(outerZmax) * m_outerTanAlpha;
  m_biggerDiscBounds = std::make_shared<RadialBounds>(
      innerRmax(), outerRmax(), averagePhi, halfPhiSector);
  m_outerConeBounds = std::make_shared<ConeBounds>(
      outerAlpha, averagePhi, outerZmin, outerZmax, halfPhiSector);
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
    innerConeRelTrans.pretranslate(
        Vector3D(0., 0., -halflengthZ() - innerOffsetZ()));
    auto innerConeAbsTrans =
        std::make_shared<Transform3D>(transform * innerConeRelTrans);
    auto innerCone =
        Surface::makeShared<ConeSurface>(innerConeAbsTrans, m_innerConeBounds);
    rSurfaces.push_back(innerCone);
  }

  // Create an outer Cone
  auto outerConeRelTrans = Transform3D::Identity();
  outerConeRelTrans.pretranslate(
      Vector3D(0., 0., -halflengthZ() - outerOffsetZ()));
  auto outerConeAbsTrans =
      std::make_shared<Transform3D>(transform * outerConeRelTrans);
  auto outerCone =
      Surface::makeShared<ConeSurface>(outerConeAbsTrans, m_outerConeBounds);
  rSurfaces.push_back(outerCone);

  // Set a disc at Zmin
  if (m_smallerDiscBounds != nullptr) {
    auto smallerDiscRelTrans = Transform3D::Identity();
    smallerDiscRelTrans.pretranslate(Vector3D(0., 0., -halflengthZ()));
    auto smallerDiscAbsTrans =
        std::make_shared<Transform3D>(transform * smallerDiscRelTrans);
    auto smallerDisc = Surface::makeShared<DiscSurface>(smallerDiscAbsTrans,
                                                        m_smallerDiscBounds);
    rSurfaces.push_back(smallerDisc);
  }

  // Set a disc at Zmax
  auto biggerDiscRelTrans = Transform3D::Identity();
  biggerDiscRelTrans.pretranslate(Vector3D(0., 0., halflengthZ()));
  auto biggerDiscAbsTrans =
      std::make_shared<Transform3D>(transform * biggerDiscRelTrans);
  auto biggerDisc =
      Surface::makeShared<DiscSurface>(biggerDiscAbsTrans, m_biggerDiscBounds);
  rSurfaces.push_back(biggerDisc);

  return rSurfaces;
}

bool Acts::ConicalVolumeBounds::inside(const Vector3D& pos, double tol) const {
  return true;
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
