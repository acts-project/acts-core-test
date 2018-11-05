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
  };

  /// @brief This struct stores the data for the construction of a cuboid
  /// TrackingVolume with a given number of PlaneLayers
  struct VolumeConfig
  {
    // Center position
    Vector3D position;
    // Lengths in x,y,z
    Vector3D length;
    // Rotation
    RotationMatrix3D rotation = RotationMatrix3D::Identity();
    // Configurations of its layers
    std::vector<LayerConfig> layerCfg;
    // Stored layers
    std::vector<std::shared_ptr<const Layer>> layers;
    // Binning direction (x,y,z)
    BinningValue binningValue;
    // Name of the volume
    std::string name = "Volume";
  };

  struct Config
  {
    std::vector<VolumeConfig>    volumeCfg;
    std::vector<TrackingVolume*> volumes;
  };

  BoxGeometryBuilder() = default;

  //~ template<typename DetectorElement_t>
  //~ std::shared_ptr<TrackingGeometry>
  //~ buildGeometry() const;

  //~ template<typename DetectorElement_t>
  //~ std::shared_ptr<TrackingGeometry>
  //~ buildGeometry(const std::vector<Vector3D>& pixelSurfaces, const
  // std::vector<Vector3D>& stripSurfaces, const std::pair<double, double>&
  // detectorLength) const;

  //~ template<typename DetectorElement_t>
  //~ std::shared_ptr<TrackingGeometry>
  //~ buildGeometry(const std::vector<Vector3D>& pixelSurfaces, const
  // std::vector<Vector3D>& stripSurfaces, const double distStrips, const
  // std::pair<double, double>& detectorLength) const;

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
  cfg.surface = buildSurface<DetectorElement_t>(cfg.surfaceCfg);
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
// TODO: treat case without surfaces
// TODO: dense volume
template <typename DetectorElement_t>
std::shared_ptr<TrackingVolume>
BoxGeometryBuilder::buildVolume(VolumeConfig& cfg) const
{
  // Build transformation
  Transform3D trafo(Transform3D::Identity() * cfg.rotation);
  trafo.translation() = cfg.position;
  // Set bounds
  auto bounds = std::make_shared<const CuboidVolumeBounds>(
      cfg.length.x() * 0.5, cfg.length.y() * 0.5, cfg.length.z() * 0.5);

  // Gather the layers
  cfg.layers.reserve(cfg.layerCfg.size());  // TODO: maybe can be removed
  LayerVector layVec;
  for (auto& layerCfg : cfg.layerCfg) {
    cfg.layers.push_back(buildLayer(layerCfg));
    layVec.push_back(cfg.layers.back());
  }

  // Build layer array
  std::pair<double, double> minMax = binningRange(cfg);
  LayerArrayCreator layArrCreator(
      getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
  std::unique_ptr<const LayerArray> layArr(
      layArrCreator.layerArray(layVec,
                               minMax.first,
                               minMax.second,
                               BinningType::arbitrary,
                               cfg.binningValue));

  // Build TrackingVolume
  auto trackVolume
      = TrackingVolume::create(std::make_shared<const Transform3D>(trafo),
                               bounds,
                               nullptr,
                               std::move(layArr),
                               layVec,
                               {},
                               {},
                               "Volume 1");
  trackVolume->sign(GeometrySignature::Global);

  return trackVolume;
}

std::pair<double, double>
BoxGeometryBuilder::binningRange(const VolumeConfig& cfg) const
{
  // Construct return value
  std::pair<double, double> minMax = std::make_pair(
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
  // Test for binning direction
  switch (cfg.binningValue) {
  // Binning in x-direction
  case BinningValue::binX: {
    for (const auto& layercfg : cfg.layerCfg) {
      // Test if new extreme is found and set it
      if (layercfg.surfaceCfg.position.x() - layercfg.layerThickness
          < minMax.first)
        minMax.first
            = layercfg.surfaceCfg.position.x() - layercfg.layerThickness;
      if (layercfg.surfaceCfg.position.x() + layercfg.layerThickness
          > minMax.second)
        minMax.second
            = layercfg.surfaceCfg.position.x() + layercfg.layerThickness;
    }
    break;
  }
  // Binning in x-direction
  case BinningValue::binY: {
    for (const auto& layercfg : cfg.layerCfg) {
      // Test if new extreme is found and set it
      if (layercfg.surfaceCfg.position.y() - layercfg.layerThickness
          < minMax.first)
        minMax.first
            = layercfg.surfaceCfg.position.y() - layercfg.layerThickness;
      if (layercfg.surfaceCfg.position.y() + layercfg.layerThickness
          > minMax.second)
        minMax.second
            = layercfg.surfaceCfg.position.y() + layercfg.layerThickness;
    }
    break;
  }
  // Binning in x-direction
  case BinningValue::binZ: {
    for (const auto& layercfg : cfg.layerCfg) {
      // Test if new extreme is found and set it
      if (layercfg.surfaceCfg.position.z() - layercfg.layerThickness
          < minMax.first)
        minMax.first
            = layercfg.surfaceCfg.position.z() - layercfg.layerThickness;
      if (layercfg.surfaceCfg.position.z() + layercfg.layerThickness
          > minMax.second)
        minMax.second
            = layercfg.surfaceCfg.position.z() + layercfg.layerThickness;
    }
    break;
  }
  // Don't do anything if direction is something else
  default: {
  }
  }
  return minMax;
}

//~ /// @brief Builds a simple 4-layer detector with 2 pixel-like and 2
//~ /// double-strip-like detectors
//~ ///
//~ /// @return Pointer to the tracking geometry
//~ template<typename DetectorElement_t>
//~ std::shared_ptr<TrackingGeometry>
//~ BoxGeometryBuilder::buildGeometry() const
//~ {
// Set translation vectors
//~ double                eps = 1. * units::_mm;
//~ std::vector<Vector3D> translations;
//~ translations.push_back({-2. * units::_m, 0., 0.});
//~ translations.push_back({-1. * units::_m, 0., 0.});
//~ translations.push_back({1. * units::_m - eps, 0., 0.});
//~ translations.push_back({1. * units::_m + eps, 0., 0.});
//~ translations.push_back({2. * units::_m - eps, 0., 0.});
//~ translations.push_back({2. * units::_m + eps, 0., 0.});

//~ // Construct surfaces
//~ std::array<PlaneSurface*, 6> surfaces;
//~ unsigned int i;
//~ for (i = 0; i < translations.size(); i++) {
//~ Transform3D trafo(Transform3D::Identity() * rotation);
//~ trafo.translation() = translations[i];

//~ surfaces[i] = new PlaneSurface(
//~ rBounds,
//~ *(new DetectorElement_t(std::make_shared<const Transform3D>(trafo),
//~ rBounds,
//~ 1. * units::_um)));
//~ surfaces[i]->setAssociatedMaterial(surMat);
//~ }

// Construct layers
//~ std::array<LayerPtr, 6> layers;
//~ for (i = 0; i < 6; i++) {
//~ Transform3D trafo(Transform3D::Identity() * rotation);
//~ trafo.translation() = translations[i];

//~ std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfaces[i]));

//~ layers[i] = PlaneLayer::create(std::make_shared<const Transform3D>(trafo),
//~ rBounds,
//~ std::move(surArray),
//~ 1. * units::_mm);
//~ surfaces[i]->associateLayer(*layers[i]);
//~ }

//~ // Build volume for surfaces with negative x-values
//~ Transform3D trafoVol1(Transform3D::Identity());
//~ trafoVol1.translation() = Vector3D(-1.5 * units::_m, 0., 0.);

//~ auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
//~ 1.5 * units::_m, 0.5 * units::_m, 0.5 * units::_m);

//~ LayerArrayCreator layArrCreator(
//~ getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
//~ LayerVector layVec;
//~ layVec.push_back(layers[0]);
//~ layVec.push_back(layers[1]);
//~ std::unique_ptr<const LayerArray> layArr1(
//~ layArrCreator.layerArray(layVec,
//~ -2. * units::_m - 1. * units::_mm,
//~ -1. * units::_m + 1. * units::_mm,
//~ BinningType::arbitrary,
//~ BinningValue::binX));

//~ auto trackVolume1
//~ = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol1),
//~ boundsVol,
//~ nullptr,
//~ std::move(layArr1),
//~ layVec,
//~ {},
//~ {},
//~ "Volume 1");
//~ trackVolume1->sign(GeometrySignature::Global);

//~ // Build volume for surfaces with positive x-values
//~ Transform3D trafoVol2(Transform3D::Identity());
//~ trafoVol2.translation() = Vector3D(1.5 * units::_m, 0., 0.);

//~ layVec.clear();
//~ for (i = 2; i < 6; i++) layVec.push_back(layers[i]);
//~ std::unique_ptr<const LayerArray> layArr2(
//~ layArrCreator.layerArray(layVec,
//~ 1. * units::_m - 2. * units::_mm,
//~ 2. * units::_m + 2. * units::_mm,
//~ BinningType::arbitrary,
//~ BinningValue::binX));

//~ auto trackVolume2
//~ = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVol2),
//~ boundsVol,
//~ nullptr,
//~ std::move(layArr2),
//~ layVec,
//~ {},
//~ {},
//~ "Volume 2");
//~ trackVolume2->sign(GeometrySignature::Global);

//~ // Glue volumes
//~ trackVolume2->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
//~ trackVolume1,
//~ BoundarySurfaceFace::positiveFaceYZ);

//~ trackVolume1->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
//~ trackVolume2,
//~ BoundarySurfaceFace::negativeFaceYZ);

//~ // Build world volume
//~ Transform3D trafoWorld(Transform3D::Identity());
//~ trafoWorld.translation() = Vector3D(0., 0., 0.);

//~ auto worldVol = std::make_shared<const CuboidVolumeBounds>(
//~ 3. * units::_m, 0.5 * units::_m, 0.5 * units::_m);

//~ std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

//~ tapVec.push_back(
//~ std::make_pair(trackVolume1, Vector3D(-1.5 * units::_m, 0., 0.)));
//~ tapVec.push_back(
//~ std::make_pair(trackVolume2, Vector3D(1.5 * units::_m, 0., 0.)));

//~ std::vector<double> binBoundaries = {-3. * units::_m, 0., 3. * units::_m};

//~ BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
//~ std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

//~ std::shared_ptr<const TrackingVolumeArray> trVolArr(
//~ new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

//~ MutableTrackingVolumePtr mtvpWorld(
//~ TrackingVolume::create(std::make_shared<const Transform3D>(trafoWorld),
//~ worldVol,
//~ trVolArr,
//~ "World"));

//~ mtvpWorld->sign(GeometrySignature::Global);

//~ // Build and return tracking geometry
//~ return std::shared_ptr<TrackingGeometry>(
//~ new Acts::TrackingGeometry(mtvpWorld));
//~ }

//~ /// @brief Builds a simple 4-layer detector with 2 pixel-like and 2
//~ /// double-strip-like detectors
//~ ///
//~ /// @return Pointer to the tracking geometry
//~ template<typename DetectorElement_t>
//~ std::shared_ptr<TrackingGeometry>
//~ BoxGeometryBuilder::buildGeometry(const std::vector<Vector3D>&
// pixelSurfaces, const std::vector<Vector3D>& stripSurfaces, const
// std::pair<double, double>& detectorLength) const
//~ {
//~ // Construct surfaces
//~ std::vector<PlaneSurface*> surfacesP =
// buildSurfaces<DetectorElement_t>(pixelSurfaces);
//~ std::vector<PlaneSurface*> surfacesS =
// buildSurfaces<DetectorElement_t>(stripSurfaces);

//~ // Construct layers
//~ std::vector<LayerPtr> layersP, layersS;
//~ layersP.reserve(pixelSurfaces.size());
//~ layersS.reserve(stripSurfaces.size());
//~ for (unsigned int i = 0; i < pixelSurfaces.size(); i++) {
//~ Transform3D trafo(Transform3D::Identity() * rotation);
//~ trafo.translation() = pixelSurfaces[i];

//~ std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfacesP[i]));

//~ layersP.push_back(PlaneLayer::create(std::make_shared<const
// Transform3D>(trafo),
//~ rBounds,
//~ std::move(surArray),
//~ 1. * units::_mm));
//~ surfacesP[i]->associateLayer(*layersP[i]);
//~ }
//~ for (unsigned int i = 0; i < stripSurfaces.size(); i++) {
//~ Transform3D trafo(Transform3D::Identity() * rotation);
//~ trafo.translation() = stripSurfaces[i];

//~ std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfacesS[i]));

//~ layersS.push_back(PlaneLayer::create(std::make_shared<const
// Transform3D>(trafo),
//~ rBounds,
//~ std::move(surArray),
//~ 1. * units::_mm));
//~ surfacesS[i]->associateLayer(*layersS[i]);
//~ }

//~ // Build volume for surfaces with negative x-values
//~ Transform3D trafoVolP(Transform3D::Identity());
//~ trafoVolP.translation() = Vector3D(detectorLength.first * 0.5, 0., 0.);

//~ auto boundsVol = std::make_shared<const CuboidVolumeBounds>(
//~ detectorLength.first, 0.5 * units::_m, 0.5 * units::_m);

//~ LayerArrayCreator layArrCreator(
//~ getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
//~ LayerVector layVec;
//~ for(unsigned int i = 0; i < layersP.size(); i++)
//~ layVec.push_back(layersP[i]);

//~ std::unique_ptr<const LayerArray> layArrP(
//~ layArrCreator.layerArray(layVec,
//~ pixelSurfaces[0].x() - 1. * units::_mm,
//~ pixelSurfaces.back().x() + 1. * units::_mm,
//~ BinningType::arbitrary,
//~ BinningValue::binX));

//~ auto trackVolumeP
//~ = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVolP),
//~ boundsVol,
//~ nullptr,
//~ std::move(layArrP),
//~ layVec,
//~ {},
//~ {},
//~ "Volume 1");
//~ trackVolumeP->sign(GeometrySignature::Global);

//~ // Build volume for surfaces with positive x-values
//~ Transform3D trafoVolS(Transform3D::Identity());
//~ trafoVolS.translation() = Vector3D(detectorLength.second * 0.5, 0., 0.);

//~ layVec.clear();
//~ for(unsigned int i = 0; i < layersS.size(); i++)
//~ layVec.push_back(layersS[i]);

//~ std::unique_ptr<const LayerArray> layArrS(
//~ layArrCreator.layerArray(layVec,
//~ stripSurfaces[0].x() - 1. * units::_mm,
//~ stripSurfaces.back().x() + 1. * units::_mm,
//~ BinningType::arbitrary,
//~ BinningValue::binX));

//~ auto trackVolumeS
//~ = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVolS),
//~ boundsVol,
//~ nullptr,
//~ std::move(layArrS),
//~ layVec,
//~ {},
//~ {},
//~ "Volume 2");
//~ trackVolumeS->sign(GeometrySignature::Global);

//~ // Glue volumes
//~ trackVolumeS->glueTrackingVolume(BoundarySurfaceFace::negativeFaceYZ,
//~ trackVolumeP,
//~ BoundarySurfaceFace::positiveFaceYZ);

//~ trackVolumeP->glueTrackingVolume(BoundarySurfaceFace::positiveFaceYZ,
//~ trackVolumeS,
//~ BoundarySurfaceFace::negativeFaceYZ);

//~ // Build world volume
//~ Transform3D trafoWorld(Transform3D::Identity());
//~ trafoWorld.translation() = Vector3D(0., 0., 0.);

//~ auto worldVol = std::make_shared<const CuboidVolumeBounds>(
//~ (detectorLength.first + detectorLength.second) * 0.5, 0.5 * units::_m, 0.5 *
// units::_m);

//~ std::vector<std::pair<TrackingVolumePtr, Vector3D>> tapVec;

//~ tapVec.push_back(
//~ std::make_pair(trackVolumeP, Vector3D(detectorLength.first * 0.5, 0., 0.)));
//~ tapVec.push_back(
//~ std::make_pair(trackVolumeS, Vector3D(detectorLength.second * 0.5, 0.,
// 0.)));

//~ std::vector<double> binBoundaries = {detectorLength.first, 0.,
// detectorLength.second};

//~ BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
//~ std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

//~ std::shared_ptr<const TrackingVolumeArray> trVolArr(
//~ new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

//~ MutableTrackingVolumePtr mtvpWorld(
//~ TrackingVolume::create(std::make_shared<const Transform3D>(trafoWorld),
//~ worldVol,
//~ trVolArr,
//~ "World"));

//~ mtvpWorld->sign(GeometrySignature::Global);

//~ // Build and return tracking geometry
//~ return std::shared_ptr<TrackingGeometry>(
//~ new Acts::TrackingGeometry(mtvpWorld));
//~ }

}  // namespace Acts
