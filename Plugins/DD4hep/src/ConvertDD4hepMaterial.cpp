// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"

void Acts::addProtoMaterial(const ActsExtension& actsExtension, Layer& layer,
                            const std::string& openBinning,
                            Acts::BinningValue openBinVal) {
  // Start with the representing surface
  std::vector<std::string> materialOptions = {"layer_material_representing"};
  std::vector<const Surface*> materialSurfaces = {
      &(layer.surfaceRepresentation())};
  // Now fill (optionally) with the approach surfaces
  auto aDescriptor = layer.approachDescriptor();
  if (aDescriptor and aDescriptor->containedSurfaces().size() >= 2) {
    // Add the inner and outer approach surface
    const std::vector<const Surface*>& aSurfaces =
        aDescriptor->containedSurfaces();
    materialOptions.push_back("layer_material_inner");
    materialSurfaces.push_back(aSurfaces[0]);
    materialOptions.push_back("layer_material_outer");
    materialSurfaces.push_back(aSurfaces[1]);
  }

  // Now loop over it and create the ProtoMaterial
  for (unsigned int is = 0; is < materialOptions.size(); ++is) {
    if (actsExtension.hasValue(materialOptions[is])) {
      unsigned int binsPhi =
          actsExtension.getValue("binsPhi", materialOptions[is]);
      unsigned int binsOpen =
          actsExtension.getValue(openBinning, materialOptions[is]);
      // Create the bin utility
      Acts::BinUtility bu;
      if (binsPhi > 1) {
        bu +=
            Acts::BinUtility(binsPhi, -M_PI, M_PI, Acts::closed, Acts::binPhi);
      }
      if (binsOpen > 1) {
        bu += Acts::BinUtility(binsOpen, -1., 1., Acts::open, openBinVal);
      }
      auto psMaterial = std::make_shared<Acts::ProtoSurfaceMaterial>(bu);
      // const_cast (ugly - to be changed after internal geometry stored
      // non-conast)
      Surface* surface = const_cast<Surface*>(materialSurfaces[is]);
      surface->assignSurfaceMaterial(psMaterial);
    }
  }
}

void Acts::addCylinderProtoMaterial(dd4hep::DetElement detElement,
                                    Layer& cylinderLayer,
                                    Logging::Level loggingLevel) {
  auto DD4hepConverterlogger =
      Acts::getDefaultLogger("DD4hepConversion", loggingLevel);
  ACTS_LOCAL_LOGGER(DD4hepConverterlogger);

  ACTS_INFO(
      "Translating DD4hep material into Acts material for CylinderLayer : "
      << detElement.name());
  // Get the Acts extension in order to prepare the ProtoMaterial
  auto actsExtension = detElement.extension<ActsExtension>();
  if (actsExtension != nullptr and actsExtension->hasType("layer_material")) {
    addProtoMaterial(*actsExtension, cylinderLayer, "binsZ", Acts::binZ);
  }
}

void Acts::xml2LayerProtoMaterial(
    const xml_comp_t& x_layer, ActsExtension& actsExtension,
    const std::vector<std::string>& materialOptions,
    const std::pair<std::string, std::string>& binOptions) {
  // Only continue if the layer has a material tag
  if (x_layer.hasChild(_Unicode(layer_material))) {
    xml_comp_t x_layer_material = x_layer.child(_Unicode(layer_material));
    // Add the layer material flag
    std::string baseTag = "layer_material";
    actsExtension.addType(baseTag);
    for (auto& materialOpt : materialOptions) {
      if (x_layer_material.attr<bool>(materialOpt.c_str())) {
        std::string materialTag = baseTag + std::string("_") + materialOpt;
        std::string bin0 = binOptions.first;
        actsExtension.addValue(x_layer_material.attr<int>(bin0.c_str()), bin0,
                               materialTag);
        std::string bin1 = binOptions.second;
        actsExtension.addValue(x_layer_material.attr<int>(bin1.c_str()), bin1,
                               materialTag);
      }
    }
  }
}

void Acts::xml2CylinderProtoMaterial(const xml_comp_t& x_layer,
                                     ActsExtension& actsExtension) {
  xml2LayerProtoMaterial(
      x_layer, actsExtension, {"inner", "representing", "outer"},
      std::pair<std::string, std::string>{"binsPhi", "binsZ"});
}

void Acts::addDiscProtoMaterial(dd4hep::DetElement detElement, Layer& discLayer,
                                Logging::Level loggingLevel) {
  auto DD4hepConverterlogger =
      Acts::getDefaultLogger("DD4hepConversion", loggingLevel);
  ACTS_LOCAL_LOGGER(DD4hepConverterlogger);

  ACTS_INFO("Translating DD4hep material into Acts material for DiscLayer : "
            << detElement.name());

  // Get the Acts extension
  auto actsExtension = detElement.extension<ActsExtension>();
  if (actsExtension != nullptr and actsExtension->hasType("layer_material")) {
    addProtoMaterial(*actsExtension, discLayer, "binsR", Acts::binR);
  }
}

void Acts::xml2DiscProtoMaterial(const xml_comp_t& x_layer,
                                 ActsExtension& actsExtension) {
  xml2LayerProtoMaterial(
      x_layer, actsExtension, {"inner", "representing", "outer"},
      std::pair<std::string, std::string>{"binsPhi", "binsR"});
}
