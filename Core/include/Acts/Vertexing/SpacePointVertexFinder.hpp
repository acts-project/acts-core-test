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
#include <utility>
#include <vector>
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;
using namespace boost::histogram;

namespace Acts {

/// @brief A projective vertex finder from space points
///
/// @tparam scanner_t is the underlying engine which performs
/// the scan and fills a target histogram

template <typename scanner_t>
class SpacePointVertexFinder {
 public:
  /// Nested configuration struct
  struct Config {
    /// The target histogram range
    std::pair<double, double> targetRange = {};
    /// The target histogram bins
    unsigned int targetBins = 300;
  };

  /// Constructor with config object and logger
  ///
  /// @param cfg is the configuration object
  /// @param logger logging instance
  SpacePointVertexFinder(const Config& cfg, scanner_t&& finder,
                         std::unique_ptr<const Logger> logger =
                             getDefaultLogger("SpacePointVertexFinder",
                                              Logging::INFO))
      : m_cfg(cfg), m_scanner(std::move(finder)), m_logger(std::move(logger)) {}

  /// Destructor
  ~SpacePointVertexFinder() = default;

  /// @brief the actual scan method
  ///
  /// @tparam seed_container_t a templated nested container
  /// @tparam seed_graph_t a templated graph descriptor
  /// @tparam confirmation_container_t a templated nested container
  ///
  /// It takes a seed container and a graph prescription and combines
  /// the seeds accordingly, then it checks with the confirmation histogram
  ///
  /// @param seedContainer the seed container for seeds building
  /// @param seedCombinations the pairwise seed combination
  /// @param confirmContainer the confirmation layer hits
  ///
  /// @return a scan histogram
  template <typename seed_container_t, typename seed_graph_t,
            typename confirmation_container_t>
  auto find(const seed_container_t& seedContainer,
            const seed_graph_t& seedGraph,
            const confirmation_container_t& confirmContainer) const {
    // Get the ordered scan lists
    auto scanLists = m_scanner.scanLists(seedContainer);

    // And get the confirmation histogram & confirmation lists ]
    auto [confirmHist, confirmLists, confirmTargets] =
        m_scanner.fillConfirmation(confirmContainer);

    // Create the target histogram with config specifications
    auto tHistogram = make_histogram(
        axis::regular<>(m_cfg.targetBins, m_cfg.targetRange.first,
                        m_cfg.targetRange.second, scanner_t::targetName()));

    // Now scan for vertices and return
    auto vertices = m_scanner.scan(tHistogram, scanLists, seedGraph,
                                   confirmHist, confirmTargets);

    // Return the target histogram wiht the vertex values
    return vertices;
  }

 private:
  /// Configuration struct
  Config m_cfg;

  /// The underlying engine
  scanner_t m_scanner;

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private acces method to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

//#include "SpacePointVertexFinder.ipp"
