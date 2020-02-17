// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

// Helper
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

// The class to test
#include "Acts/Geometry/Extent.hpp"

namespace Acts {

using namespace UnitLiterals;

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(ExtentTest) {

    std::vector<Vector3D> vertices = 
        { Vector3D(15_mm,-3_mm, -10_mm),
          Vector3D(18_mm, 0_mm, -10_mm),
          Vector3D(15_mm,-3_mm, -10_mm),
          Vector3D(15_mm,-3_mm,  10_mm),
          Vector3D(18_mm, 0_mm,  10_mm),
          Vector3D(15_mm,-3_mm,  10_mm)};

    // Create an Extent
    Extent gExt;
    for (const auto& v : vertices){
        gExt.check(v);
    }

    CHECK_CLOSE_ABS(gExt.ranges[binX].first,15_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binX].second,18_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binY].first,-3_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binY].second,3_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binZ].first,0_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binZ].second,0_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binR].first,15_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binR].second,18_mm,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binPhi].first, ,1e-6);
    CHECK_CLOSE_ABS(gExt.ranges[binPhi].second, ,1e-6);

}

BOOST_AUTO_TEST_SUITE_END()

}
}