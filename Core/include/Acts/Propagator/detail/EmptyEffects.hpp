// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace detail {

constexpr int number_of_component = 2;

struct EmptyEffect {
  struct ComponentValues {
    double weight;
    double mean;
    double variance;
  };
  std::vector<ComponentValues> getMixture(double /*unused*/,
                                          double /*unused*/) const {
    // make non effect, just means copy the components in material effect
    std::vector<ComponentValues> comp(number_of_component,
                                      {1. / number_of_component, 0., 0.});
    return std::move(comp);
  }
};
}  // namespace detail
}  // namespace Acts
