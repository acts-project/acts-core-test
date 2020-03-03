// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

/// @class ProtoSurfaceMaterial
///
/// @brief proxy to SurfaceMaterial hand over BinUtility
///
/// The ProtoSurfaceMaterial class acts as a proxy to the SurfaceMaterial
/// to mark the layers and surfaces on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity of
/// of the material map with the bin Utility.

class ProtoVolumeMaterial : public IVolumeMaterial {
 public:
  /// Constructor without BinUtility - homogenous material
  ProtoVolumeMaterial() = default;

  /// Constructor with BinUtility - multidimensional material
  ///
  /// @param binUtility a BinUtility determining the granularity
  ///        and binning of the material on the surface/layer
  ProtoVolumeMaterial(const BinUtility& binUtility);

  /// Copy constuctor
  ///
  /// @param smproxy The source proxy
  ProtoVolumeMaterial(const ProtoVolumeMaterial& smproxy) = default;

  /// Copy move constuctor
  ///
  /// @param smproxy The source proxy
  ProtoVolumeMaterial(ProtoVolumeMaterial&& smproxy) = default;

  /// Destructor
  ///
  /// @param smproxy The source proxy
  ~ProtoVolumeMaterial() override = default;

  /// Assignment operator
  ///
  /// @param smproxy The source proxy
  ProtoVolumeMaterial& operator=(const ProtoVolumeMaterial& smproxy) = default;

  /// Return the material
  const Material& material(const Vector3D& /*position*/) const final {
    return m_material;
  }

 private:
  Material m_material;
};
}  // namespace Acts
