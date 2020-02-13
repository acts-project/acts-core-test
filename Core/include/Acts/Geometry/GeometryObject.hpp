// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

using range_type = std::pair<double, double>;

/// @class GeometryObject
///
/// Base class to provide GeometryID interface:
/// - simple set and get
///
/// It also provides the binningPosition method for
/// Geometry geometrical object to be binned in BinnedArrays
///
class GeometryObject {
 public:
  /// @brief Extent in space
  ///
  /// This is a nested struct to the GeometryObject representation
  /// which can be retrieved and used for surface parsing and will
  /// give you the maximal extent in 3D space/
  struct Extent {
    /// Possible maximal value
    static constexpr double maxval = std::numeric_limits<double>::max();

    /// The range in x
    range_type xrange = {maxval, -maxval};
    /// The range in x
    range_type yrange = {maxval, -maxval};
    /// The range in x
    range_type zrange = {maxval, -maxval};
    /// The range in x
    range_type rrange = {maxval, -maxval};

    /// Check the vertex
    /// @param vtx the Vertex to be checked
    void check(const Vector3D& vtx) {
      // min/max value check
      auto minMax = [&](range_type& range, double value) -> void {
        range.first = std::min(value, range.first);
        range.second = std::max(value, range.second);
      };
      // go through the parameters
      minMax(xrange, vtx.x());
      minMax(yrange, vtx.y());
      minMax(zrange, vtx.z());
      minMax(rrange, VectorHelpers::perp(vtx));
    }
  };

  /// Defaulted construrctor
  GeometryObject() = default;

  /// Defaulted copy constructor
  GeometryObject(const GeometryObject&) = default;

  /// Constructor from a ready-made value
  ///
  /// @param geoID the geometry identifier of the object
  GeometryObject(const GeometryID& geoID) : m_geoID(geoID) {}

  /// Assignment operator
  ///
  /// @param geoID the source geoID
  GeometryObject& operator=(const GeometryObject& geoID) {
    if (&geoID != this) {
      m_geoID = geoID.m_geoID;
    }
    return *this;
  }

  /// @return the geometry id by reference
  const GeometryID& geoID() const;

  /// Force a binning position method
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the value in which you want to bin
  ///
  /// @return vector 3D used for the binning schema
  virtual const Vector3D binningPosition(const GeometryContext& gctx,
                                         BinningValue bValue) const = 0;

  /// Implement the binningValue
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the dobule in which you want to bin
  ///
  /// @return float to be used for the binning schema
  virtual double binningPositionValue(const GeometryContext& gctx,
                                      BinningValue bValue) const;

  /// Set the value
  ///
  /// @param geoID the geometry identifier to be assigned
  void assignGeoID(const GeometryID& geoID);

 protected:
  GeometryID m_geoID;
};

inline const GeometryID& GeometryObject::geoID() const {
  return m_geoID;
}

inline void GeometryObject::assignGeoID(const GeometryID& geoID) {
  m_geoID = geoID;
}

inline double GeometryObject::binningPositionValue(const GeometryContext& gctx,
                                                   BinningValue bValue) const {
  using VectorHelpers::perp;
  // now switch
  switch (bValue) {
    // case x
    case Acts::binX: {
      return binningPosition(gctx, bValue).x();
    } break;
    // case y
    case Acts::binY: {
      return binningPosition(gctx, bValue).y();
    } break;
    // case z
    case Acts::binZ: {
      return binningPosition(gctx, bValue).z();
    } break;
    // case - to be overwritten by disks
    case Acts::binR: {
      return perp(binningPosition(gctx, bValue));
    } break;
    // do nothing for the default
    default:
      return 0.;
  }
}
}  // namespace Acts
