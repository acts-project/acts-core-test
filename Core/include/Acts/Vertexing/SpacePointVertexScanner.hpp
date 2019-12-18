// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <boost/histogram.hpp>
#include <cmath>
#include <tuple>
#include <vector>
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/SpacePointVertexProjectors.hpp"

using namespace Acts::UnitLiterals;
using namespace boost::histogram;

namespace Acts {

using Bin = std::pair<unsigned int, unsigned int>;

/// This is the input class
template <typename input_hit_t>
struct ScanPoint {
  /// Link to the input object
  input_hit_t input;

  /// Prepped data : rho
  double rho = 0.;

  /// Prepped data : phi
  double phi = 0.;

  /// The test bin assignment
  unsigned int itb = 0;

  // The hit bin raw values
  Eigen::Array<double, 1, 3> rawValues;

  /// The hit bin assignments
  Bin bin = {0, 0};

  /// Data point where you copy
  ScanPoint(const input_hit_t& idata) : input(idata) {}

  /// Data point where you move it into
  ScanPoint(input_hit_t&& idata) : input(std::move(idata)) {}
};

template <unsigned int OBINS, typename projecter_t, typename confirm_targets_t>
class SpacePointVertexScanner {
 public:
  /// Nested Config struct
  struct Config {
    /// How many do you require for acceptance
    unsigned int require = 1;
    /// Binning
    unsigned int confirmBins0 = 0;
    unsigned int confirmBins1 = 0;
    /// The vertex target
    typename confirm_targets_t::value_type vertexTarget;
    /// Cell mask - for additional bin setting
    std::vector<Bin> neighborMask = {};
  };

  /// Constructor with config object and logger
  ///
  /// @param cfg is the configuration object
  /// @param logger logging instance
  SpacePointVertexScanner(const Config& cfg,
                          std::unique_ptr<const Logger> logger =
                              getDefaultLogger("SpacePointVertexScanner",
                                               Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  /// Create the ordered Lists for the processing
  ///
  /// @tparam seed_container_t the templated input container
  ///
  /// @param seedContainer is a nested container of seed spacepoints
  ///
  /// @return ordered lists
  template <typename seed_container_t>
  auto scanLists(const seed_container_t& seedContainer) const {
    // Deduce and declare the type
    using InputVector = typename seed_container_t::value_type;
    using InputData = typename InputVector::value_type;
    using BinVector = std::vector<ScanPoint<InputData>>;
    using ScanList = std::array<BinVector, OBINS>;

    // Instantiate the scan lists and fill them with prepared data
    std::vector<ScanList> scanLists;
    for (const auto& seedVector : seedContainer) {
      // The ordered list for this container entry
      ScanList olist;
      for (const auto& seed : seedVector) {
        // Create the scan point for this seed & prepare
        ScanPoint<InputData> spoint(std::move(seed));
        spoint.phi = VectorHelpers::phi(seed.position());
        spoint.itb = m_tbins.index(spoint.phi);
        olist[spoint.itb].push_back(spoint);
      }
      scanLists.push_back(olist);
    }
    return scanLists;
  }

  /// Create the confirmation histogram
  ///
  /// @return a boost histogram
  template <typename confirmation_container_t>
  auto fillConfirmation(const confirmation_container_t& confirmations) const {
    // The confirmation Histogram
    auto confirmHist =
        make_histogram(axis::circular<>(m_cfg.confirmBins0, -M_PI, M_PI, "phi"),
                       axis::regular<>(m_cfg.confirmBins1, -1., 1., "z|rho"));

    // Deduce and declare the type
    using InputVector = typename confirmation_container_t::value_type;
    using InputData = typename InputVector::value_type;
    using ConfirmCollection = std::vector<ScanPoint<InputData>>;
    using ConfirmContainer = std::vector<ConfirmCollection>;

    // That's the type of the confirmation target
    using ConfirmationTarget = typename confirm_targets_t::value_type;

    // Min / max values
    double vmin = std::numeric_limits<double>::min();
    double vmax = std::numeric_limits<double>::max();

    // One for raw value casting - uses static method
    ConfirmationTarget ctarget;
    std::vector<unsigned int> ctargetIDs;

    // Preparation loop for the confirmation hits
    ConfirmContainer confirmContainer;
    confirmContainer.reserve(confirmations.size());
    // Prepration of the min/max values for the targets
    std::vector<Eigen::Array<double, 2, 2>> mmMatrix;
    mmMatrix.reserve(confirmations.size());
    // Loop and parse
    for (const auto& confirmRow : confirmations) {
      ConfirmCollection confirmCollection;
      confirmCollection.reserve(confirmRow.size());
      // The min/max array
      Eigen::Array<double, 2, 2> mm;
      mm << vmin, vmax, vmin, vmax;
      // Fill the raw values and get min max
      for (const auto& confirmHit : confirmRow) {
        ScanPoint<InputData> spoint(std::move(confirmHit));
        spoint.rawValues = ctarget.rawValues(spoint.input.position());
        mm(0, 0) = std::min(mm(0, 0), spoint.rawValues(0, 1));
        mm(0, 1) = std::max(mm(0, 1), spoint.rawValues(0, 1));
        mm(1, 0) = std::min(mm(1, 0), spoint.rawValues(0, 2));
        mm(1, 1) = std::max(mm(1, 1), spoint.rawValues(0, 2));
        confirmCollection.push_back(std::move(spoint));
      }
      mmMatrix.push_back(mm);
      confirmContainer.push_back(confirmCollection);
    }

    // Create the confirmation targets per container & rescale
    confirm_targets_t confirmTargets;
    size_t ic = 0;
    for (auto& confirmCollection : confirmContainer) {
      // Get the min/max parsing
      confirmTargets.push_back(ConfirmationTarget(mmMatrix[ic], ++ic));
      const auto& target = confirmTargets.back();
      for (auto& confirmHit : confirmCollection) {
        auto pvals = target.projectionParameters(confirmHit);
        target.scaleProjection(pvals);
        confirmHit.bin = {confirmHist.axis(0).index(pvals(0, 0)),
                          confirmHist.axis(1).index(pvals(0, 1))};

        // Get the value, assign the hit and put it back
        std::vector<Bin> binMask = m_cfg.neighborMask;
        binMask.push_back(Bin{0, 0});
        unsigned int ib0 = confirmHit.bin.first;
        unsigned int ib1 = confirmHit.bin.second;
        for (const auto& bmsk : binMask) {
          int im0 = ib0 + bmsk.first;
          int im1 = ib1 + bmsk.second;
          // Under / overlow protection for regular binning
          if (im1 < 0 or im1 > m_cfg.confirmBins1) {
            continue;
          }
          // Wrap around for circular binning
          im0 = im0 < 0 ? m_cfg.confirmBins0 - 1
                        : im0 + 1 > m_cfg.confirmBins0 ? 0 : im0;
          // Now get and set
          uint64_t bval = confirmHist.at(im0, im1);
          bval |= (1 << target.targetID());
          confirmHist.at(im0, im1) = bval;
        }
      }
    }

    // Return the histogram and the massaged data
    return std::tuple<decltype(confirmHist), ConfirmContainer,
                      confirm_targets_t>{std::move(confirmHist),
                                         std::move(confirmContainer),
                                         std::move(confirmTargets)};
  }

  template <typename vertex_hist_t, typename seed_lists_t,
            typename seed_graph_t, typename confirm_hist_t>
  auto scan(vertex_hist_t& vertices, const seed_lists_t& seeds,
            const seed_graph_t& seedGraph, const confirm_hist_t& confirmHist,
            const confirm_targets_t& confirmTargets) const {
    // Run over the seed combination / graph if you wanna name it so
    for (auto sgraph : seedGraph) {
      // The first is a reference list
      auto refIter = sgraph.begin();
      auto nextIter = refIter;
      const auto& refList = seeds[*refIter];
      while (++nextIter != sgraph.end()) {
        const auto& nextList = seeds[*nextIter];
        auto refBinIter = refList.begin();
        // Iterator access to get the synchronised list
        for (; refBinIter < refList.end(); ++refBinIter) {
          auto ibin = std::distance(refList.begin(), refBinIter);
          // Get the synchronised container
          const auto& nextBin = nextList[ibin];
          for (const auto& ref : (*refBinIter)) {
            for (const auto& next : nextBin) {
              // Look for confirmations
              unsigned int acceptCount = 0;
              auto projz = m_projector.project(ref, next, m_cfg.vertexTarget);
              const auto& axis0 = confirmHist.axis(0);
              const auto& axis1 = confirmHist.axis(1);
              // if (projz(0,1) < or projz(0,1) > ) continue;
              for (const auto& ctarget : confirmTargets) {
                auto proj = m_projector.project(ref, next, ctarget);
                auto ib0 = axis0.index(double(proj(0, 0)));
                auto ib1 = axis1.index(double(proj(0, 1)));
                uint64_t bval = confirmHist.at(ib0, ib1);
                // Check if the bit is set for this value
                if (bval & (1 << (ctarget.targetID()))) {
                  ++acceptCount > m_cfg.require;
                  vertices(projz(0, 1));
                }
              }
            }
          }
        }
      }
    }
    return vertices;
  }

  /// The target name
  static const std::string targetName() { return std::string("z"); }

 private:
  /// The config struct
  Config m_cfg;

  /// The projector
  projecter_t m_projector;

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private acces method to the logger
  const Logger& logger() const { return *m_logger; }

  /// The test axis for the engine
  axis::circular<> m_tbins = axis::circular<>{OBINS, -M_PI, M_PI, "tphi"};
};

}  // namespace Acts
