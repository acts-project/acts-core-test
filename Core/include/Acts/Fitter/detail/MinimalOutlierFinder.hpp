// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/math/distributions/chi_squared.hpp>
#include "Acts/EventData/TrackState.hpp"

namespace Acts {

struct MinimalOutlierFinder {
  /// The measurement significance criteria
  std::map<OutlierSearchStage, double> measurementSignificance;

  /// @brief Public call mimicking an outlier searcher
  ///
  /// @param surface The surface of the measurement
  /// @param chi2 The chisq from fitting
  /// @param ndf The measurement dimension
  /// @param searchStage The outlier search stage
  ///
  /// @return The resulting
  bool operator()(const Surface* /*surface*/, double chi2, size_t ndf,
                  OutlierSearchStage searchStage) const {
    if (measurementSignificance.find(searchStage) !=
        measurementSignificance.end()) {
      // The chisq distribution
      boost::math::chi_squared chiDist(ndf);
      // The p-Value
      double pValue = 1 - boost::math::cdf(chiDist, chi2);
      // If pValue is NOT significant enough => outlier
      return pValue > measurementSignificance.at(searchStage) ? false : true;
    }
    return false;
  }
};

}  // namespace Acts
