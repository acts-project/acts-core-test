// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class SurfaceArray;
using CompoundSurfaceArrays = BinnedArray<std::unique_ptr<SurfaceArray>>;
using CompoundApproachDescriptors =
    BinnedArray<std::unique_ptr<ApproachDescriptor>>;

// @class CompoundLayer
///
/// Class to describe a Layer that is compound by sub layers
/// for tracking, it still has a single surface representation which
/// is primarily existing for ordering the layer into a layer array
/// of  a volume.
///
class CompoundLayer : public Layer {
 public:
  /// Factory for shared Layer pointer
  ///
  /// @param repSurface The representing surface
  /// @param repThickness The representing thickness
  /// @param surfaceArrays The sub surface arrays
  /// @param approachDescriptors The individual approach descriptors
  /// @param laytyp The layer type indicator
  ///
  /// @return The return object is a shared poiter to the layer.
  static MutableLayerPtr create(
      std::shared_ptr<Surface> repSurface, double repThickness,
      std::unique_ptr<CompoundSurfaceArrays> surfaceArrays,
      std::unique_ptr<CompoundApproachDescriptors> apporachDescriptors,
      LayerType laytyp = active) {
    return MutableLayerPtr(new CompoundLayer(
        std::move(repSurface), repThickness, std::move(surfaceArrays),
        std::move(apporachDescriptors), laytyp));
  }

  /// Copy constructor - deleted
  CompoundLayer(const CompoundLayer& cla) = delete;

  /// Assignment operator for CylinderLayers - deleted
  CompoundLayer& operator=(const CompoundLayer&) = delete;

  /// Default Constructor - deleted
  CompoundLayer() = delete;

  /// Destructor - defaulted
  ~CompoundLayer() override = default;

  /// Transforms the layer into a Surface representation
  /// This is for positioning and extrapolation
  const Surface& surfaceRepresentation() const final;

  // Non-const version
  Surface& surfaceRepresentation() final;

  /// The binning position method - this will call the binning position
  /// method of the representing surface
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the type of global binning to be done
  ///
  /// @return is the global position to be used for binning
  const Vector3D binningPosition(const GeometryContext& gctx,
                                 BinningValue bValue) const final;

  /// Geometrical isOnLayer() method
  ///
  /// @note for a compound layer is checking to be within either of the
  /// representing volumes
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the gobal position to be checked
  /// @param bcheck is the boundary check directive
  ///
  /// @return boolean that indicates success of the operation
  bool isOnLayer(const GeometryContext& gctx, const Vector3D& position,
                 const BoundaryCheck& bcheck = true) const final;

  /// Return the corresponding SurfaceArray, returns a nullptr if no
  /// SurfaceArray
  ///
  /// @param position is an optional argument if multiple surface arrays
  /// exist in case of a compound layer
  const SurfaceArray* surfaceArray(const Vector3D& position) const final;

  /// Return method for the approach descriptor, can be nullptr
  ///
  /// @param position is an optional argument if multiple
  /// approach descriptors exist, e.g. for compount layers
  const ApproachDescriptor* approachDescriptor(
      const Vector3D& position) const final;

  /// Return a Polyhedron for this object
  ///
  /// Will return a compound polyhedron built from the sub components
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, if the lseg
  /// is set to one, only the corners and the extrema are given,
  /// otherwise it represents the number of segments for a full 2*M_PI
  /// circle and is scaled to the relevant sector
  ///
  /// @note An internal surface transform can invalidate the extrema
  /// in the transformed space
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const final;

 private:
  std::shared_ptr<Surface> m_representingSurface;
  std::unique_ptr<CompoundSurfaceArrays> m_surfaceArrays;
  std::unique_ptr<CompoundApproachDescriptors> m_apporachDescriptors;

 protected:
  /// Private constructor for CompoundLayer
  ///
  /// @param repSurface The representing surface
  /// @param repThickness The representing thickness
  /// @param surfaceArrays The sub surface arrays
  /// @param approachDescriptors The individual approach descriptors
  /// @param laytyp The layer type indicator
  ///
  /// @return The return object is a shared poiter to the layer.
  CompoundLayer(
      std::shared_ptr<Surface> repSurface, double repThickness,
      std::unique_ptr<CompoundSurfaceArrays> surfaceArrays,
      std::unique_ptr<CompoundApproachDescriptors> apporachDescriptors,
      LayerType laytyp = active);
};

inline const SurfaceArray* CompoundLayer::surfaceArray(
    const Vector3D& position) const {
  if (m_surfaceArrays != nullptr) {
    std::array<size_t, 3> bin;
    m_surfaceArrays->object(position, bin).get();
  }
  return nullptr;
}

inline const ApproachDescriptor* CompoundLayer::approachDescriptor(
    const Vector3D& position) const {
  if (m_apporachDescriptors != nullptr) {
    std::array<size_t, 3> bin;
    m_apporachDescriptors->object(position, bin).get();
  }
  return nullptr;
}

inline const Surface& CompoundLayer::surfaceRepresentation() const {
  return (*(m_representingSurface.get()));
}

inline Surface& CompoundLayer::surfaceRepresentation() {
  return (*(m_representingSurface.get()));
}

bool CompoundLayer::isOnLayer(const GeometryContext& gctx,
                              const Vector3D& position,
                              const BoundaryCheck& bcheck) const {
  return true;
}

inline const Vector3D CompoundLayer::binningPosition(
    const GeometryContext& gctx, BinningValue bValue) const {
  return m_representingSurface->binningPosition(gctx, bValue);
}

}  // namespace Acts