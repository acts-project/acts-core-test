// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;
class CylinderBounds;
class ConeBounds;
class RadialBounds;
class PlanarBounds;
class Volume;

/// @class ConicalVolumeBounds
///
/// Bounds for a cubical Volume, the decomposeToSurfaces method creates a
/// vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///

class ConicalVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for readability
  enum BoundValues {
    bv_innerAlpha = 0,
    bv_innerOffsetZ = 1,
    bv_outerAlpha = 2,
    bv_outerOffsetZ = 3,
    bv_halflengthZ = 4,
    bv_averagePhi = 5,
    bv_halfPhiSector = 6,
    bv_length = 7
  };

  /// Deleted default Constructor
  ConicalVolumeBounds() = delete;

  /// Constructor - for general cone-cone setups
  ///
  /// @param innerAlpha The opening angle of the inner cone (0 if no cone)
  /// @param innerOffsetZ The tip  z position in of the inner cone, w.r.t center
  /// @param outerAlpha  The opening angle of the outer cone (0 if no cone)
  /// @param outerOffsetZ The tip  z position in of the outer cone, w.r.t center
  /// @param halflengthZ The minimum z value of the inner and outer cones
  /// @param avergePhi The phi orientation of the sector
  /// @param halfPhiSector The opening angle phi sector
  ConicalVolumeBounds(double innerAlpha, double innerTipZ, double outerAlpha,
                      double outerOffsetZ, double halflengthZ,
                      double averagePhi, double halfPhiSector);

  /// Constructor - for general cylidner-cone setups
  ///
  /// @param cylinderR The inner radius of the cylinder
  /// @param alpha  The opening angle of the cone (0 if no cone)
  /// @param offsetZ The tip  z position in of the cone, w.r.t center
  /// @param halflengthZ The minimum z value of the inner and outer cones
  /// @param avergePhi The phi orientation of the sector (defaulted to 0)
  /// @param halfPhiSector The opening angle phi sector
  ///
  /// @note depending on cylinderR > coneR it is constructing a cone with
  /// cylindrical cutout or a cylinder with conical cutout
  ConicalVolumeBounds(double cylinderR, double alpha, double offsetZ,
                      double halflengthZ, double averagePhi,
                      double halfPhiSector);

  /// Copy Constructor
  ///
  /// @param cobo is the source volume bounds to be copied
  ConicalVolumeBounds(const ConicalVolumeBounds& cobo);

  /// Defaulted destructor
  ~ConicalVolumeBounds() override = default;

  /// Assignment operator
  ///
  /// @param cobo is the source volume bounds to be assigned
  ConicalVolumeBounds& operator=(const ConicalVolumeBounds& cobo);

  /// Virtual constructor
  ConicalVolumeBounds* clone() const override;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  ///
  /// @param pos is the position in volume frame to be checked
  /// @param tol is the absolute tolerance to be applied
  bool inside(const Vector3D& pos, double tol = 0.) const final;

  /// Method to decompose the Bounds into boundarySurfaces
  ///
  /// @param transformPtr is the transfrom of the volume
  SurfacePtrVector decomposeToSurfaces(
      const Transform3D* transformPtr = nullptr) const final;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3D* trf = nullptr,
                                  const Vector3D& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Return the inner alpha value
  VOLUMEBOUNDS_VALUESTORE_ACCESS(innerAlpha);

  /// Return the z value of the inner tip
  VOLUMEBOUNDS_VALUESTORE_ACCESS(innerOffsetZ);

  /// Return the outer alpha value
  VOLUMEBOUNDS_VALUESTORE_ACCESS(outerAlpha);

  /// Return the z value of the outer tip
  VOLUMEBOUNDS_VALUESTORE_ACCESS(outerOffsetZ);

  /// Return the halflength in Z
  VOLUMEBOUNDS_VALUESTORE_ACCESS(halflengthZ);

  /// Return the average phi value
  VOLUMEBOUNDS_VALUESTORE_ACCESS(averagePhi);

  /// Return the half phi sector value
  VOLUMEBOUNDS_VALUESTORE_ACCESS(halfPhiSector);

  // Return the derived inner tan(alpha)
  VOLUMEBOUNDS_DERIVED_ACCESS(innerTanAlpha);

  // Return the derived innerRmin
  VOLUMEBOUNDS_DERIVED_ACCESS(innerRmin);

  // Return the derived innerRmin
  VOLUMEBOUNDS_DERIVED_ACCESS(innerRmax);

  // Return the derived inner tan(alpha)
  VOLUMEBOUNDS_DERIVED_ACCESS(outerTanAlpha);

  // Return the derived innerRmin
  VOLUMEBOUNDS_DERIVED_ACCESS(outerRmin);

  // Return the derived innerRmin
  VOLUMEBOUNDS_DERIVED_ACCESS(outerRmax);

  /// Output Method for std::ostream
  ///
  /// @param sl is ostream operator to be dumped into
  std::ostream& toStream(std::ostream& sl) const override;

 private:
  // Create the sector bounds if needed
  void createSectorBounds();

  /// Templated dumpT method
  template <class T>
  T& dumpT(T& dt) const;

  /// The bound values
  std::vector<TDD_real_t> m_valueStore;
  std::shared_ptr<CylinderBounds> m_innerCylinderBounds = nullptr;
  std::shared_ptr<ConeBounds> m_innerConeBounds = nullptr;
  std::shared_ptr<ConeBounds> m_outerConeBounds = nullptr;
  std::shared_ptr<CylinderBounds> m_outerCylinderBounds = nullptr;
  std::shared_ptr<RadialBounds> m_negativeDiscBounds = nullptr;
  std::shared_ptr<RadialBounds> m_positiveDiscBounds = nullptr;
  std::shared_ptr<PlanarBounds> m_sectorBounds = nullptr;

  /// Derived values
  double m_innerRmin = 0.;
  double m_innerRmax = 0.;
  double m_outerRmin = 0.;
  double m_outerRmax = 0.;

  double m_innerTanAlpha = 0.;
  double m_outerTanAlpha = 0.;
};

template <class T>
T& ConicalVolumeBounds::dumpT(T& dt) const {
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(5);
  dt << "Acts::ConicalVolumeBounds :";

  return dt;
}
}  // namespace Acts