// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>
#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"

namespace Acts {

/// @brief This class builds a box detector with a configurable amount of
/// surfaces in it.
class BoxGeometryBuilder
{
public:
  /// @brief This struct stores the data for the construction of a single
  /// PlaneSurface
  struct SurfaceConfig
  {
    // Center position
    Vector3D position;
    // Rotation
    RotationMatrix3D rotation = RotationMatrix3D::Identity();
    // Bounds
    std::shared_ptr<const RectangleBounds> rBounds = nullptr;
    // Attached material
    std::shared_ptr<const SurfaceMaterial> surMat = nullptr;
    // Thickness
    double thickness = 0.;
  };

  /// @brief This struct stores the data for the construction of a PlaneLayer
  /// that has a single PlaneSurface encapsulated
  struct LayerConfig
  {
    // Configuration of the surface
    SurfaceConfig surfaceCfg;
    // Thickness
    double layerThickness = 0.;
    // Encapsulated surface
    PlaneSurface* surface = nullptr;
    // Boolean flag if layer is active
    bool active = false;
  };

  /// @brief This struct stores the data for the construction of a cuboid
  /// TrackingVolume with a given number of PlaneLayers
  struct VolumeConfig
  {
    // Center position
    Vector3D position;
    // Lengths in x,y,z
    Vector3D length;
    // Configurations of its layers
    std::vector<LayerConfig> layerCfg;
    // Stored layers
    std::vector<std::shared_ptr<const Layer>> layers;
    // Name of the volume
    std::string name = "Volume";
    // Material
    std::shared_ptr<const Material> material = nullptr;
  };

  /// @brief This struct stores the configuration of the tracking geometry
  struct Config
  {
    // Center position
    Vector3D position;
    // Length in x,y,z
    Vector3D length;
    // Configuration of its volumes
    std::vector<VolumeConfig> volumeCfg;
    // Stored volumes
    std::vector<std::shared_ptr<TrackingVolume>> volumes;
  };

  BoxGeometryBuilder() = default;

  /// @brief This function creates a surface with a given configuration. A
  /// detector element is attached if the template parameter is non-void.
  ///
  /// @tparam DetectorElement_t Type of the optional detector element
  /// @param [in] cfg Configuration of the surface
  /// @return Pointer to the created surface
  template <typename DetectorElement_t = void>
  PlaneSurface*
  buildSurface(const SurfaceConfig& cfg) const;

  /// @brief This function creates a layer with a surface encaspulated with a
  /// given configuration. The surface gets a detector element attached if the
  /// template parameter is non-void.
  ///
  /// @tparam DetectorElement_t Type of the optional detector element
  /// @param [in, out] cfg Configuration of the layer and the surface
  /// @return Pointer to the created layer
  template <typename DetectorElement_t = void>
  std::shared_ptr<const Layer>
  buildLayer(LayerConfig& cfg) const;

  /// @brief This function creates a TrackingVolume with a configurable number
  /// of layers and surfaces. Each surface gets a detector element attached if
  /// the template parameter is non-void.
  ///
  /// @tparam DetectorElement_t Type of the optional detector element
  /// @param [in, out] cfg Configuration of the TrackingVolume
  /// @return Pointer to the created TrackingVolume
  template <typename DetectorElement_t = void>
  std::shared_ptr<TrackingVolume>
  buildVolume(VolumeConfig& cfg) const;

  /// @brief This function builds a TrackingGeometry based on a given
  /// configuration
  ///
  /// @tparam DetectorElement_t Type of the optional detector element
  /// @param [in, out] cfg Configuration of the TrackingGeometry
  /// @return Pointer to the created TrackingGeometry
  template <typename DetectorElement_t = void>
  std::shared_ptr<TrackingGeometry>
  buildTrackingGeometry(Config& cfg) const;

  /// @brief This function evaluates the minimum and maximum of the binning as
  /// given by the configurations of the surfaces and layers. The ordering
  /// depends on the binning value specified in the configuration of the volume.
  ///
  /// @param [in] cfg Container with the given surfaces and layers
  /// @return Pair containing the minimum and maximum along the binning
  /// direction
  std::pair<double, double>
  binningRange(const VolumeConfig& cfg) const;
};

template <typename DetectorElement_t>
PlaneSurface*
BoxGeometryBuilder::buildSurface(const SurfaceConfig& cfg) const
{
  PlaneSurface* surface;

  // Build transformation
  Transform3D trafo(Transform3D::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;

  // Create and store surface
  surface = new PlaneSurface(
      cfg.rBounds,
      *(new DetectorElement_t(std::make_shared<const Transform3D>(trafo),
                              cfg.rBounds,
                              cfg.thickness)));
  surface->setAssociatedMaterial(cfg.surMat);
  return surface;
}

template <>
Acts::PlaneSurface*
Acts::BoxGeometryBuilder::buildSurface<void>(const SurfaceConfig& cfg) const
{
  PlaneSurface* surface;

  // Build transformation
  Transform3D trafo(Transform3D::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;

  // Create and store surface
  surface = new PlaneSurface(std::make_shared<const Transform3D>(trafo),
                             cfg.rBounds);
  surface->setAssociatedMaterial(cfg.surMat);
  return surface;
}

template <typename DetectorElement_t>
std::shared_ptr<const Layer>
BoxGeometryBuilder::buildLayer(LayerConfig& cfg) const
{
  LayerPtr layer;

  // Build the surface
  if (cfg.surface == nullptr) {
    cfg.surface = buildSurface<DetectorElement_t>(cfg.surfaceCfg);
  }
  // Build transformation centered at the surface position
  Transform3D trafo(Transform3D::Identity() * cfg.surfaceCfg.rotation);
  trafo.translation() = cfg.surfaceCfg.position;

  // Get the surface and build the layer
  std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(cfg.surface));

  layer
      = PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
                           cfg.surfaceCfg.rBounds,
                           std::move(surArray),
                           cfg.layerThickness,
                           nullptr,
                           (cfg.surface->associatedDetectorElement() == nullptr
                                ? LayerType::passive
                                : LayerType::active));
  cfg.surface->associateLayer(*layer);
  return layer;
}

template <typename DetectorElement_t>
std::shared_ptr<TrackingVolume>
BoxGeometryBuilder::buildVolume(VolumeConfig& cfg) const
{
  // Build transformation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = cfg.position;
  // Set bounds
  auto bounds = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  if (cfg.layerCfg.empty()) {
    // Build dummy layer if no layer is given (tmp solution)
    SurfaceConfig sCfg;
    sCfg.position = cfg.position;
    // Rotation of the surfaces
    double   rotationAngle = M_PI * 0.5;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    sCfg.rotation.col(0) = xPos;
    sCfg.rotation.col(1) = yPos;
    sCfg.rotation.col(2) = zPos;
    // Bounds
    sCfg.rBounds = std::make_shared<const RectangleBounds>(
        RectangleBounds(cfg.length.y() * 0.5, cfg.length.z() * 0.5));

    LayerConfig lCfg;
    lCfg.surfaceCfg     = sCfg;
    lCfg.layerThickness = 1. * units::_mm;

    cfg.layerCfg.push_back(lCfg);
  }

  // Gather the layers
  LayerVector layVec;
  if (cfg.layers.empty()) {
    cfg.layers.reserve(cfg.layerCfg.size());

    for (auto& layerCfg : cfg.layerCfg) {
      if (layerCfg.active) {
        cfg.layers.push_back(buildLayer<DetectorElement_t>(layerCfg));
      } else {
        cfg.layers.push_back(buildLayer(layerCfg));
      }
      layVec.push_back(cfg.layers.back());
    }
  } else {
    for (auto& lay : cfg.layers) {
      layVec.push_back(lay);
    }
  }

  // Build layer array
  std::pair<double, double> minMax = binningRange(cfg);
  LayerArrayCreator layArrCreator(
      getDefaultLogger("LayerArrayCreator", Logging::INFO));
  std::unique_ptr<const LayerArray> layArr(
      layArrCreator.layerArray(layVec,
                               minMax.first,
                               minMax.second,
                               BinningType::arbitrary,
                               BinningValue::binX));

  // Build TrackingVolume
  auto trackVolume
      = TrackingVolume::create(std::make_shared<const Transform3D>(trafo),
                               bounds,
                               cfg.material,
                               std::move(layArr),
                               layVec,
                               {},
                               {},
                               cfg.name);
  trackVolume->sign(GeometrySignature::Global);

  return trackVolume;
}

std::pair<double, double>
BoxGeometryBuilder::binningRange(const VolumeConfig& cfg) const
{
  // Construct return value
  std::pair<double, double> minMax = std::make_pair(
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
  for (const auto& layercfg : cfg.layerCfg) {
    // Test if new extreme is found and set it
    if (layercfg.surfaceCfg.position.x() - layercfg.layerThickness
        < minMax.first) {
      minMax.first = layercfg.surfaceCfg.position.x() - layercfg.layerThickness;
    }
    if (layercfg.surfaceCfg.position.x() + layercfg.layerThickness
        > minMax.second) {
      minMax.second
          = layercfg.surfaceCfg.position.x() + layercfg.layerThickness;
    }
  }
  return minMax;
}

template <typename DetectorElement_t>
std::shared_ptr<TrackingGeometry>
BoxGeometryBuilder::buildTrackingGeometry(Config& cfg) const
{
  // Build volumes
  if (cfg.volumes.empty()) {
    cfg.volumes.reserve(cfg.volumeCfg.size());
    for (VolumeConfig volCfg : cfg.volumeCfg) {
      cfg.volumes.push_back(buildVolume<DetectorElement_t>(volCfg));
    }
  }

  // Glue volumes
  // TODO: YZ due to x-binning. Keep it that way or allow variations?
  for (unsigned int i = 0; i < cfg.volumes.size() - 1; i++) {
    cfg.volumes[i + 1]->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
                                           cfg.volumes[i],
                                           BoundarySurfaceFace::positiveFaceYZ);
    cfg.volumes[i]->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
                                       cfg.volumes[i + 1],
                                       BoundarySurfaceFace::negativeFaceYZ);
  }

  // Translation
  Transform3D trafo(Transform3D::Identity());
  trafo.translation() = cfg.position;

  // Size of the volume
  auto volume = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  // Build vector of confined volumes
  std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;
  tapVec.reserve(cfg.volumeCfg.size());
  for (auto& tVol : cfg.volumes) {
    tapVec.push_back(std::make_pair(tVol, tVol->center()));
  }

  // Set bin boundaries along binning
  std::vector<float> binBoundaries;
  binBoundaries.push_back(cfg.volumes[0]->center().x()
                          - cfg.volumeCfg[0].length.x() * 0.5);

  for (size_t i = 0; i < cfg.volumes.size(); i++) {
    binBoundaries.push_back(cfg.volumes[i]->center().x()
                            + cfg.volumeCfg[i].length.x() * 0.5);
  }

  // Build binning
  BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
  std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

  // Build TrackingVolume array
  std::shared_ptr<const TrackingVolumeArray> trVolArr(
      new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

  // Create world volume
  MutableTrackingVolumePtr mtvp(TrackingVolume::create(
      std::make_shared<const Transform3D>(trafo), volume, trVolArr, "World"));

  mtvp->sign(GeometrySignature::Global);

  // Build and return tracking geometry
  return std::shared_ptr<TrackingGeometry>(new TrackingGeometry(mtvp));
}
}  // namespace Acts
