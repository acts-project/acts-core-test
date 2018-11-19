// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE CuboidVolumeBuilderTest

#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

  ///
  /// @brief Stub implementation of the detector element
  ///
  class DetElem : public DetectorElementBase
  {
  public:
    /// @brief Constructor
    ///
    /// @param [in] trafo Transformation of the detector element
    /// @param [in] rBounds Rectangle boundaries of the plane surface
    /// @param [in] thickness Thickness of the detector element
    DetElem(std::shared_ptr<const Transform3D>     trafo,
            std::shared_ptr<const RectangleBounds> rBounds,
            double                                 thickness)
      : DetectorElementBase()
      , m_trafo(trafo)
      , m_surface(new PlaneSurface(rBounds, *this))
      , m_thickness(thickness)
    {
    }

    /// @brief Getter of the transformation
    virtual const Transform3D&
    transform() const
    {
      return *m_trafo;
    }

    /// @brief Getter of the surface
    virtual const Surface&
    surface() const
    {
      return *m_surface;
    }

    /// @brief Getter of the thickness
    virtual double
    thickness() const
    {
      return m_thickness;
    }

    // Pointer to the transformation
    std::shared_ptr<const Transform3D> m_trafo;
    // Surface related to the detector element
    Surface const* m_surface;
    // Thickness of the detector element
    double m_thickness;
  };

  BOOST_AUTO_TEST_CASE(CuboidVolumeBuilderTest)
  {

    // Construct builder
    CuboidVolumeBuilder cvb;

    // Create configurations for surfaces
    std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
    for (unsigned int i = 1; i < 5; i++) {
      // Position of the surfaces
      CuboidVolumeBuilder::SurfaceConfig cfg;
      cfg.position = {i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;

      surfaceConfig.push_back(cfg);
    }

    // Test that there are actually 4 surface configurations
    BOOST_TEST(surfaceConfig.size() == 4);

    // Test that 4 sensitive surfaces can be built
    for (const auto& cfg : surfaceConfig) {
      PlaneSurface* pSur = cvb.buildSurface<DetElem>(cfg);
      BOOST_TEST(pSur != nullptr);
      BOOST_TEST(pSur->center() == cfg.position);
      BOOST_TEST(pSur->associatedMaterial() != nullptr);
      BOOST_TEST(pSur->associatedDetectorElement() != nullptr);
    }

    // Test that 4 passive surfaces can be built
    for (const auto& cfg : surfaceConfig) {
      PlaneSurface* pSur = cvb.buildSurface<>(cfg);
      BOOST_TEST(pSur != nullptr);
      BOOST_TEST(pSur->center() == cfg.position);
      BOOST_TEST(pSur->associatedMaterial() != nullptr);
      BOOST_TEST(pSur->associatedDetectorElement() == nullptr);
    }

    ////////////////////////////////////////////////////////////////////
    // Build layer configurations
    std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig;
    for (auto& sCfg : surfaceConfig) {
      CuboidVolumeBuilder::LayerConfig cfg;
      cfg.surfaceCfg = sCfg;
      layerConfig.push_back(cfg);
    }

    // Test that there are actually 4 layer configurations
    BOOST_TEST(layerConfig.size() == 4);

    // Test that 4 layers with sensitive surfaces can be built
    for (auto& cfg : layerConfig) {
      LayerPtr layer = cvb.buildLayer<DetElem>(cfg);
      BOOST_TEST(layer != nullptr);
      BOOST_TEST(cfg.surface != nullptr);
      BOOST_TEST(layer->surfaceArray()->surfaces().size() == 1);
      BOOST_TEST(layer->layerType() == LayerType::active);
    }

    // Test that 4 layers with passive surfaces can be built
    for (auto& cfg : layerConfig) {
      cfg.surface    = nullptr;
      LayerPtr layer = cvb.buildLayer<>(cfg);
      BOOST_TEST(layer != nullptr);
      BOOST_TEST(cfg.surface != nullptr);
      BOOST_TEST(layer->surfaceArray()->surfaces().size() == 1);
    }
    for (auto& cfg : layerConfig) {
      cfg.surface = nullptr;
    }

    // Build volume configuration
    CuboidVolumeBuilder::VolumeConfig volumeConfig;
    volumeConfig.position = {2.5 * units::_m, 0., 0.};
    volumeConfig.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig.layerCfg = layerConfig;
    volumeConfig.name     = "Test volume";
    volumeConfig.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));

    // Test the building
    std::shared_ptr<TrackingVolume> trVol = cvb.buildVolume(volumeConfig);
    BOOST_TEST(volumeConfig.layers.size() == 4);
    BOOST_TEST(trVol->confinedLayers()->arrayObjects().size()
               == volumeConfig.layers.size() * 2
                   + 1);  // #layers = navigation + material layers
    BOOST_TEST(trVol->volumeName() == volumeConfig.name);
    BOOST_TEST(trVol->material() != nullptr);

    // Test the building
    volumeConfig.layers.clear();
    trVol = cvb.buildVolume<DetElem>(volumeConfig);
    BOOST_TEST(volumeConfig.layers.size() == 4);
    BOOST_TEST(trVol->confinedLayers()->arrayObjects().size()
               == volumeConfig.layers.size() * 2
                   + 1);  // #layers = navigation + material layers
    BOOST_TEST(trVol->volumeName() == volumeConfig.name);

    volumeConfig.layers.clear();
    for (auto& lay : volumeConfig.layerCfg) {
      lay.surface = nullptr;
      lay.active  = true;
    }
    trVol = cvb.buildVolume<DetElem>(volumeConfig);
    BOOST_TEST(volumeConfig.layers.size() == 4);
    for (auto& lay : volumeConfig.layers) {
      BOOST_TEST(lay->layerType() == LayerType::active);
    }
    
    volumeConfig.layers.clear();
    for (auto& lay : volumeConfig.layerCfg) {
      lay.active = true;
    }
    trVol = cvb.buildVolume<DetElem>(volumeConfig);
    BOOST_TEST(volumeConfig.layers.size() == 4);
    for (auto& lay : volumeConfig.layers) {
      BOOST_TEST(lay->layerType() == LayerType::active);
    }

    ////////////////////////////////////////////////////////////////////
    // Build TrackingGeometry configuration

    // Build second volume
    std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig2;
    for (int i = 1; i < 5; i++) {
      // Position of the surfaces
      CuboidVolumeBuilder::SurfaceConfig cfg;
      cfg.position = {-i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;
      surfaceConfig2.push_back(cfg);
    }

    std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig2;
    for (auto& sCfg : surfaceConfig2) {
      CuboidVolumeBuilder::LayerConfig cfg;
      cfg.surfaceCfg = sCfg;
      layerConfig2.push_back(cfg);
    }
    CuboidVolumeBuilder::VolumeConfig volumeConfig2;
    volumeConfig2.position = {-2.5 * units::_m, 0., 0.};
    volumeConfig2.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig2.layerCfg = layerConfig2;
    volumeConfig2.name     = "Test volume2";

    CuboidVolumeBuilder::Config config;
    config.position  = {0., 0., 0.};
    config.length    = {10. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {volumeConfig2, volumeConfig};

    cvb.setConfig(config);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder tgb(tgbCfg);

    std::unique_ptr<const TrackingGeometry> detector = tgb.trackingGeometry();
    BOOST_TEST(detector->lowestTrackingVolume({1., 0., 0.})->volumeName()
               == volumeConfig.name);
    BOOST_TEST(detector->lowestTrackingVolume({-1., 0., 0.})->volumeName()
               == volumeConfig2.name);
  }

  BOOST_AUTO_TEST_CASE(BoxGeometryBuilderTest_confinedVolumes)
  {
    // Production factory
    BoxGeometryBuilder bgb;

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg;
    vCfg.position = {1. * units::_m, 0., 0.};
    vCfg.length   = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg.name     = "Test volume";
    // Build and add 2 confined volumes
    BoxGeometryBuilder::VolumeConfig cvCfg1;
    cvCfg1.position = {1.1 * units::_m, 0., 0.};
    cvCfg1.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg1.name     = "Confined volume1";
    cvCfg1.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    BoxGeometryBuilder::VolumeConfig cvCfg2;
    cvCfg2.position = {0.9 * units::_m, 0., 0.};
    cvCfg2.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg2.name     = "Confined volume2";
    vCfg.volumeCfg  = {cvCfg1, cvCfg2};

    // Build detector
    BoxGeometryBuilder::Config config;
    config.position  = {1. * units::_m, 0., 0.};
    config.length    = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {vCfg};
    std::shared_ptr<TrackingGeometry> detector
        = bgb.buildTrackingGeometry(config);

    // Test that the right volume is selected
    BOOST_TEST(
        detector->lowestTrackingVolume({1. * units::_m, 0., 0.})->volumeName()
        == vCfg.name);
    BOOST_TEST(
        detector->lowestTrackingVolume({1.1 * units::_m, 0., 0.})->volumeName()
        == cvCfg1.name);
    BOOST_TEST(
        detector->lowestTrackingVolume({0.9 * units::_m, 0., 0.})->volumeName()
        == cvCfg2.name);

    // Set propagator and navigator
    PropagatorOptions<ActionList<StepVolumeCollector>> propOpts;
    propOpts.maxStepSize = 10. * units::_mm;
    StraightLineStepper sls;
    Navigator           navi(detector);
    navi.resolvePassive   = true;
    navi.resolveMaterial  = true;
    navi.resolveSensitive = true;

    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        nullptr, startParams, startMom, 1.);

    // Launch and collect results
    const auto& result = prop.propagate(sbtp, propOpts);
    const StepVolumeCollector::this_result& stepResult
        = result.get<typename StepVolumeCollector::result_type>();

    // Check the identified volumes
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (i > 0) {
        BOOST_TEST(stepResult.position[i].x() > 0.);
      }
      if (stepResult.position[i].x() >= 0.85 * units::_m
          && stepResult.position[i].x() < 0.95 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg2.name);
        BOOST_TEST(stepResult.volume[i]->material() == nullptr);
      } else if (stepResult.position[i].x() >= 1.05 * units::_m
                 && stepResult.position[i].x() < 1.15 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg1.name);
        BOOST_TEST(stepResult.volume[i]->material() != nullptr);
      } else if (stepResult.position[i].x() < 2. * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg.name);
        BOOST_TEST(stepResult.volume[i]->material() == nullptr);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(BoxGeometryBuilderTest_confinedVolumes_edgecases)
  {
    // Production factory
    BoxGeometryBuilder bgb;

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg1;
    vCfg1.position = {1. * units::_m, 0., 0.};
    vCfg1.length   = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg1.name     = "Test volume1";
    // Build and add 4 confined volumes
    // Volume that is missed and quite orthogonal to the starting position
    BoxGeometryBuilder::VolumeConfig cvCfg1;
    cvCfg1.position = {0.1 * units::_m, 0.4 * units::_m, 0.4 * units::_m};
    cvCfg1.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg1.name     = "Confined volume1";
    cvCfg1.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    // Volume that is missed but far away such that it may be hit
    BoxGeometryBuilder::VolumeConfig cvCfg2;
    cvCfg2.position = {1.9 * units::_m, -0.4 * units::_m, -0.4 * units::_m};
    cvCfg2.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg2.name     = "Confined volume2";
    // Volume that is hit but with identical boundary as its mother
    BoxGeometryBuilder::VolumeConfig cvCfg3;
    cvCfg3.position = {1.95 * units::_m, 0., 0.};
    cvCfg3.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg3.name     = "Confined volume3";
    // Volume to grind along the boundary
    BoxGeometryBuilder::VolumeConfig cvCfg4;
    cvCfg4.position = {1. * units::_m, 5. * units::_cm, 0.};
    cvCfg4.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg4.name     = "Confined volume4";
    vCfg1.volumeCfg = {cvCfg1, cvCfg2, cvCfg3, cvCfg4};

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg2;
    vCfg2.position = {2.5 * units::_m, 0., 0.};
    vCfg2.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg2.name     = "Test volume2";

    // Build detector
    BoxGeometryBuilder::Config config;
    config.position  = {1.5 * units::_m, 0., 0.};
    config.length    = {3. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {vCfg1, vCfg2};
    std::shared_ptr<TrackingGeometry> detector
        = bgb.buildTrackingGeometry(config);

    // Set propagator and navigator
    PropagatorOptions<ActionList<StepVolumeCollector>> propOpts;
    propOpts.maxStepSize = 10. * units::_mm;
    StraightLineStepper sls;
    Navigator           navi(detector);
    navi.resolvePassive   = true;
    navi.resolveMaterial  = true;
    navi.resolveSensitive = true;

    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        nullptr, startParams, startMom, 1.);

    // Launch and collect results
    const auto& result = prop.propagate(sbtp, propOpts);
    const StepVolumeCollector::this_result& stepResult
        = result.get<typename StepVolumeCollector::result_type>();

    // Check the identified volumes
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (i > 0) {
        BOOST_TEST(stepResult.position[i].x() > 0.);
      }
      if (stepResult.position[i].x() >= 0.95 * units::_m
          && stepResult.position[i].x() < 1.05 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg4.name);
      } else {
        if (stepResult.position[i].x() >= 1.9 * units::_m
            && stepResult.position[i].x() < 2. * units::_m) {
          BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg3.name);
        } else {
          if (stepResult.position[i].x() < 2. * units::_m) {
            BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg1.name);
          } else {
            if (stepResult.position[i].x() < 3. * units::_m)
              BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg2.name);
          }
        }
      }
    }
  }
}  // namespace Test
}  // namespace Acts
