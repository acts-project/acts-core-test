// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackState.hpp"

namespace Acts {

struct MinimalOutlierFinder {
  /// The outlier search criteria
  std::map<OutlierSearchStage, double> outlierCriteria;

  /// @brief Public call mimicking an outlier rejector
  ///
  /// @tparam chi2 The chisq from fitting
  ///
  /// @param surface The surface of the measurement
  /// @param searchStage The outlier search stage
  ///
  /// @return The resulting
  //  template <typename track_state_t>
  bool operator()(double chi2, const Surface* /*surface*/,
                  OutlierSearchStage searchStage) const {
    if (outlierCriteria.find(searchStage) != outlierCriteria.end()) {
      return chi2 > outlierCriteria.at(searchStage) ? true : false;
    }
    return false;
  }
};

}  // namespace Acts
