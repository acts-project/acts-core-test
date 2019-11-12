// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/range/adaptors.hpp>
#include <memory>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/covariance_helper.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
template <typename parameters_t>
class GainMatrixSmoother {
  using jacobian_t = typename parameters_t::CovMatrix_t;

 public:
  /// @brief Gain Matrix smoother implementation
  ///

  /// Constructor with (non-owning) logger
  /// @param logger a logger instance
  GainMatrixSmoother(
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("GainMatrixSmoother", Logging::INFO).release()))
      : m_logger(std::move(logger)) {}

  template <typename track_states_t>
  Result<parameters_t> operator()(
      const GeometryContext& gctx, track_states_t& filteredStates,
      const OutlierFinder outlierFinder = nullptr) const {
    ACTS_VERBOSE("Invoked GainMatrixSmoother");
    using namespace boost::adaptors;

    using track_state_t = typename track_states_t::value_type;
    using ParVector_t = typename parameters_t::ParVector_t;
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    using gain_matrix_t = CovMatrix_t;

    // smoothed parameter vector and covariance matrix
    ParVector_t smoothedPars;
    CovMatrix_t smoothedCov;

    // Find the smoothing starting state: smoothed is filtered - also: switch to
    // next
    ACTS_VERBOSE("Getting previous track state");
    track_state_t* prev_ts = nullptr;
    typename track_states_t::reverse_iterator lastFiltered = std::find_if(
        filteredStates.rbegin(), filteredStates.rend(), [](const auto& fst) {
          if (fst.parameter.filtered &&
              !fst.isType(TrackStateFlag::OutlierFlag))
            return true;
          return false;
        });

    if (lastFiltered == filteredStates.rend()) {
      return KalmanFitterError::SmoothFailed;
    }

    prev_ts = &(*lastFiltered);
    prev_ts->parameter.smoothed = *prev_ts->parameter.filtered;

    // Smoothing gain matrix
    gain_matrix_t G;

    size_t dist = std::distance(filteredStates.rbegin(), lastFiltered);
    // Loop and smooth the remaining states
    for (track_state_t& ts : filteredStates |
                                 sliced(0, filteredStates.size() - dist - 1) |
                                 reversed) {
      // Skip smoothing for the state if it's not filtered or an outlier
      if (!ts.parameter.filtered or ts.isType(TrackStateFlag::OutlierFlag)) {
        continue;
      }

      // The current state
      assert(ts.parameter.filtered);
      assert(ts.parameter.predicted);
      assert(ts.parameter.jacobian);
      assert(ts.parameter.predicted->covariance());
      assert(ts.parameter.filtered->covariance());

      assert(prev_ts->parameter.smoothed);
      assert(prev_ts->parameter.predicted);

      ACTS_VERBOSE("Calculate smoothing matrix:");
      ACTS_VERBOSE("Filtered covariance:\n"
                   << *ts.parameter.filtered->covariance());
      ACTS_VERBOSE("Jacobian:\n" << ts.parameter.jacobian->transpose());
      ACTS_VERBOSE("Prev. predicted covariance\n"
                   << *prev_ts->parameter.predicted->covariance()
                   << "\n, inverse: \n"
                   << (*prev_ts->parameter.predicted->covariance()).inverse());

      // clang-format off

      // Gain smoothing matrix
      G = (*ts.parameter.filtered->covariance())
          * ts.parameter.jacobian->transpose()
          * (*prev_ts->parameter.predicted->covariance()).inverse();
      ACTS_VERBOSE("Gain smoothing matrix is:\n" << G);

      // Calculate the smoothed parameters

        ACTS_VERBOSE("Calculate smoothed parameters:");
        ACTS_VERBOSE("Filtered parameters: " << ts.parameter.filtered->parameters().transpose());
        ACTS_VERBOSE(
            "Prev. smoothed parameters: " << prev_ts->parameter.smoothed->parameters().transpose());
        ACTS_VERBOSE(
            "Prev. predicted parameters: " << prev_ts->parameter.predicted->parameters().transpose());

      if (G.hasNaN()) {
        return  KalmanFitterError::SmoothFailed;
      }

      smoothedPars = ts.parameter.filtered->parameters()
                     + G * (prev_ts->parameter.smoothed->parameters()
                            - prev_ts->parameter.predicted->parameters());
      ACTS_VERBOSE("Smoothed parameters are: " << smoothedPars.transpose());

      // And the smoothed covariance

        ACTS_VERBOSE("Calculate smoothed covariance:");
        ACTS_VERBOSE("Prev. smoothed covariance:\n"
                     << *prev_ts->parameter.smoothed->covariance());

      smoothedCov = (*ts.parameter.filtered->covariance())
                    - G * (*(prev_ts->parameter.predicted->covariance())
                           - (*prev_ts->parameter.smoothed->covariance()))
                           * G.transpose();
      ACTS_VERBOSE("Smoothed covariance is: \n" << smoothedCov);
     
      // Check if the covariance matrix is semi-positive definite.
      // If not, attempt to replace it with the nearest semi-positive def matrix.
      if(not detail::covariance_helper<CovMatrix_t>::validate(smoothedCov)){
        ACTS_DEBUG("Non semi-positive smoothed covariance found. Smoothing failed for this state.");
        continue;
      }

      // clang-format on

      // Create smoothed track parameters
      ts.parameter.smoothed =
          parameters_t(gctx, smoothedCov, smoothedPars,
                       ts.referenceSurface().getSharedPtr());

      // Update the chi2 using smoothed track parameters and covariance
      bool isChi2Updated = false;
      if (ts.measurement.calibrated) {
        isChi2Updated = std::visit(
            [&](const auto& calibrated) {
              // type of measurement
              using meas_t =
                  typename std::remove_const<typename std::remove_reference<
                      decltype(calibrated)>::type>::type;
              // measurement covariance matrix
              using meas_cov_t = typename meas_t::CovMatrix_t;
              // measurement (local) parameter vector
              using meas_par_t = typename meas_t::ParVector_t;
              // type of projection
              using projection_t = typename meas_t::Projection_t;

              ACTS_VERBOSE("Measurement dimension: " << meas_t::size());
              ACTS_VERBOSE("Calibrated measurement: "
                           << calibrated.parameters().transpose());
              ACTS_VERBOSE("Calibrated measurement covariance:\n"
                           << calibrated.covariance());

              // Take the projector (measurement mapping function)
              const projection_t& H = calibrated.projector();
              ACTS_VERBOSE("Measurement projector H:\n" << H);

              // The residual from smoothed parameter
              meas_par_t residual = calibrated.residual(*ts.parameter.smoothed);
              ACTS_VERBOSE("Residual w.r.t. smoothed parameter: "
                           << residual.transpose());

              // The residual covariance inverse
              meas_cov_t resCovInv =
                  (calibrated.covariance() - H * smoothedCov * H.transpose())
                      .inverse();
              ACTS_VERBOSE("Residual covariance invserse: " << resCovInv);

              // Check if the covariance matrix is semi-positive definite.
              // If not, attempt to replace it with the nearest semi-positive
              // def matrix.
              if (not detail::covariance_helper<meas_cov_t>::validate(
                      resCovInv)) {
                ACTS_DEBUG(
                    "Residual covariance inverse is not semi-positive. Chi2 is "
                    "not updated.");
                return false;
              }

              // Overwrite the chi2
              ts.parameter.chi2 =
                  (residual.transpose() * resCovInv * residual).value();
              ACTS_VERBOSE("Chi2 after smoothing: " << ts.parameter.chi2);
              return true;
            },
            *ts.measurement.calibrated);
      }

      // If the chi2 is updated successfully with smoothed parameter,
      // proceed to check if the measurement is an outlier
      bool isOutlier = false;
      if (isChi2Updated && outlierFinder) {
        isOutlier = outlierFinder(&ts.referenceSurface(), ts.parameter.chi2,
                                  *ts.size(), OutlierSearchStage::Smoothing);
      }

      // Point prev state to current state if current state is NOT an outlier
      if (not isOutlier) {
        prev_ts = &ts;
      } else {
        ACTS_VERBOSE(
            "The state is found to be outlier during smoothing. Prev. state "
            "not updated.");
        ts.setType(TrackStateFlag::OutlierFlag, true);
        continue;
      }
    }
    // The result is the pointer to the last smoothed state - for the cache
    return *prev_ts->parameter.smoothed;
  }

  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const {
    assert(m_logger);
    return *m_logger;
  }
};
}  // namespace Acts
