// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file SolenoidBFieldTests.cpp
#define BOOST_TEST_MODULE Solenoid magnetic field tests


// clang-format off
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
// clang-format on

#include <fstream>


namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {
  
  BOOST_AUTO_TEST_CASE(TestSolenoidBField)
  {
    SolenoidBField::Config cfg;
    cfg.L = 5.8 * Acts::units::_m;
    cfg.R = (2.56+2.46)*0.5*0.5 * Acts::units::_m;
    cfg.nCoils = 1154;
    cfg.bMagCenter = 2. * Acts::units::_T;
    SolenoidBField bField(cfg);

    SolenoidBField::Cache cache;
    BOOST_TEST(bField.getField({0, 0, 0}, cache).isApprox(Vector3D(0, 0, 2.0 * Acts::units::_T)));

    //std::ofstream outf("solenoid.csv");
    //outf << "x;y;z;B_x;B_y;B_z" << std::endl;

    double tol = 1e-6;
    double tol_B = 1e-6 * Acts::units::_T;
    size_t steps = 20;
    for(size_t i=0;i<steps;i++) {
      double r = 1.5*cfg.R/steps * i;
      BOOST_TEST_CONTEXT("r=" << r) {
        //std::cout << cfg.R << " " << steps << " " << cfg.R/steps << " " << i << " " << r << std::endl;
        Vector3D B1 = bField.getField({r, 0, 0}, cache);
        Vector3D B2 = bField.getField({-r, 0, 0}, cache);
        BOOST_CHECK_SMALL(B1.x(), tol);
        BOOST_CHECK_SMALL(B1.y(), tol);
        BOOST_TEST(std::abs(B1.z()) > tol_B); // greater than zero
        // check symmetry: at z=0 it should be exactly symmetric
        BOOST_TEST(B1.isApprox(B2));

        // at this point in r, go along the length
        for(size_t j=0;j<=steps;j++) {
          //double z = cfg.L/steps * j - (cfg.L/2.);
          double z = (1.5*cfg.L/2.)/steps * j;
          BOOST_TEST_CONTEXT("z=" << z) {
            Vector3D B_zp_rp = bField.getField({r, 0, z}, cache);
            Vector3D B_zn_rp = bField.getField({r, 0, -z}, cache);
            Vector3D B_zp_rn = bField.getField({-r, 0, z}, cache);
            Vector3D B_zn_rn = bField.getField({-r, 0, -z}, cache);

            //outf << r << ";0;" << z << ";" << B_zp_rp.x() << ";" << B_zp_rp.y() << ";" << B_zp_rp.z() << std::endl;
            //if(j>0) {
              //outf << r << ";0;" << -z << ";" << B_zn_rp.x() << ";" << B_zn_rp.y() << ";" << B_zn_rp.z() << std::endl;
            //}
            //if(i>0) {
              //outf << -r << ";0;" << z << ";" << B_zp_rn.x() << ";" << B_zp_rn.y() << ";" << B_zp_rn.z() << std::endl;
            //}
            //if(i>0 && j>0) {
              //outf << -r << ";0;" << -z << ";" << B_zn_rn.x() << ";" << B_zn_rn.y() << ";" << B_zn_rn.z() << std::endl;
            //}

            // non-zero z
            BOOST_TEST(std::abs(B_zp_rp.z()) > tol_B);
            BOOST_TEST(std::abs(B_zn_rp.z()) > tol_B);
            BOOST_TEST(std::abs(B_zn_rn.z()) > tol_B);
            BOOST_TEST(std::abs(B_zp_rn.z()) > tol_B);
            if(i>0) {
              // z components should be the same for +- r
              BOOST_CHECK_CLOSE(B_zp_rp.z(), B_zp_rn.z(), tol);
              BOOST_CHECK_CLOSE(B_zn_rp.z(), B_zn_rn.z(), tol);
              // x components should be exactly opposite
              BOOST_CHECK_CLOSE(B_zp_rp.x(), -B_zp_rn.x(), tol);
              BOOST_CHECK_CLOSE(B_zn_rp.x(), -B_zn_rn.x(), tol);
            }
            if(j>0) {
              // z components should be the same for +- z
              BOOST_CHECK_CLOSE(B_zp_rp.z(), B_zn_rp.z(), tol);
              BOOST_CHECK_CLOSE(B_zp_rn.z(), B_zn_rn.z(), tol);
              // x components should be exactly opposite
              BOOST_CHECK_CLOSE(B_zp_rp.x(), -B_zn_rp.x(), tol);
              BOOST_CHECK_CLOSE(B_zp_rn.x(), -B_zn_rn.x(), tol);
            }
          }
        }
      }
    }
    //outf.close();
  }

}  // namespace Test

}  // namespace Acts
