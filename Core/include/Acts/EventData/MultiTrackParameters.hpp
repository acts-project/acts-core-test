// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <list>
#include <map>
#include <type_traits>
#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @class MultiTrackParameters
///
/// @brief base class for a multi set of track parameters
///
/// this class represents a multi set of track parameters (used in Gsf)

///
/// @tparam ChargePolicy type for distinguishing charged and neutral
/// tracks/particles (must be either ChargedPolicy or NeutralPolicy)
template <class ChargePolicy, class ParametersT>
class MultiTrackParameters : public TrackParametersBase {
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

 public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = typename TrackParametersBase::ParVector_t;

  /// type of covariance matrix
  using CovMatrix_t = typename TrackParametersBase::CovMatrix_t;

  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;

  /// The TrackMap sort the track by weight automatically
  using ParameterWeightTrack = std::pair<double, ParametersT>;
  using ParameterMapWeightTrack =
      std::multimap<double, ParametersT, std::greater<double>>;

  using multiSurfaceType = typename ParametersT::surfaceType;

  /// @brief standard constructor for track parameters of charged particles
  /// param list of BoundParameters
  template <typename P = ParametersT, typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters(std::list<std::pair<double, P>>&& weightParList,
                       std::shared_ptr<multiSurfaceType> surface)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_mSurface(std::move(surface)),
        m_oChargePolicy(weightCharge()),
        m_oParameters(weightParameterSet()),
        m_oTime(weightTime()),
        m_vPosition(weightPosition()),
        m_vMomentum(weightMomentum()) {
    assert(m_mSurface);
  }

  /// @brief standard constructor for track parameters of charged particles
  /// param list of CurvilinearParameters
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters(
      std::list<
          std::pair<double, SingleCurvilinearTrackParameters<ChargedPolicy>>>&&
          weightParList)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_oChargePolicy(weightCharge()),
        m_oParameters(weightParameterSet()),
        m_oTime(weightTime()),
        m_vPosition(weightPosition()),
        m_vMomentum(weightMomentum()) {
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param list of BoundParameters
  template <typename P = ParametersT, typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters(std::list<std::pair<double, P>>&& weightParList,
                       std::shared_ptr<multiSurfaceType> surface)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_mSurface(std::move(surface)),
        m_oChargePolicy(weightCharge()),
        m_oParameters(weightParameterSet()),
        m_oTime(weightTime()),
        m_vPosition(weightPosition()),
        m_vMomentum(weightMomentum()) {
    assert(m_mSurface);
  }
  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param list of CurvilinearParameters
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters(
      std::list<
          std::pair<double, SingleCurvilinearTrackParameters<ChargedPolicy>>>&&
          weightParList)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_oChargePolicy(weightCharge()),
        m_oParameters(weightParameterSet()),
        m_oTime(weightTime()),
        m_vPosition(weightPosition()),
        m_vMomentum(weightMomentum()) {
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// @brief copy constructor
  MultiTrackParameters(
      const MultiTrackParameters<ChargePolicy, ParametersT>& copy) = default;

  /// @brief default move constructor
  MultiTrackParameters(MultiTrackParameters<ChargePolicy, ParametersT>&& copy) =
      default;

  /// @brief assignment operator for BoundParameters
  MultiTrackParameters<ChargePolicy, ParametersT>& operator=(
      const MultiTrackParameters<ChargePolicy, ParametersT>& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = rhs.m_oChargePolicy;
      m_oTime = rhs.m_oTime;
      m_oParameters = rhs.m_oParameters;
      m_vPosition = rhs.m_vPosition;
      m_vMomentum = rhs.m_vMomentum;
      m_mWeightTracks = rhs.m_mWeightTracks;
      m_mSurface = rhs.m_mSurface;
    }
    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @param rhs object to be movied into `*this`
  MultiTrackParameters<ChargePolicy, ParametersT>& operator=(
      MultiTrackParameters<ChargePolicy, ParametersT>&& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = std::move(rhs.m_oChargePolicy);
      m_oTime = std::move(rhs.m_oTime);
      m_oParameters = std::move(rhs.m_oParameters);
      m_vPosition = std::move(rhs.m_vPosition);
      m_vMomentum = std::move(rhs.m_vMomentum);
      m_mWeightTracks = std::move(rhs.m_mWeightTracks);
      m_mSurface = std::move(rhs.m_mSurface);
    }
    return *this;
  }

  /// @brief default virtual destructor
  ~MultiTrackParameters() final = default;

  /// @brief equality operator
  ///
  bool operator==(const TrackParametersBase& /*unused*/) const final {
    return true;
  }

  /// @brief append a single parameter to the track map for BoundParameters
  template <typename P = ParametersT>
  void append(double weight, P&& para) {
    // They should be on the same surface
    assert(&para.referenceSurface() == m_mSurface.get());
    m_mWeightTracks.insert(std::make_pair(weight, std::move(para)));
    m_oChargePolicy = weightCharge();
    m_oParameters = getParameterSet();
    m_oTime = weightTime();
    m_vPosition = weightPosition();
    m_vMomentum = weightMomentum();
  }

  /// @brief append a single parameter to the track map for
  /// CurvilinearParameters
  void append(double weight,
              SingleCurvilinearTrackParameters<ChargedPolicy>&& para) {
    m_mWeightTracks.insert(std::make_pair(weight, std::move(para)));
    m_oChargePolicy = weightCharge();
    m_oParameters = getParameterSet();
    m_oTime = weightTime();
    m_vPosition = weightPosition();
    m_vMomentum = weightMomentum();
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// getter functions for position
  ActsVectorD<3> position() const final { return m_vPosition; }
  /// getter functions for momentum
  ActsVectorD<3> momentum() const final { return m_vMomentum; }
  /// getter functions for time
  double time() const final { return m_oTime; }
  /// getter functions for charge
  double charge() const final { return m_oChargePolicy.getCharge(); }

  /// getter functions for parset
  const FullParameterSet& getParameterSet() const final {
    return m_oParameters;
  }

  /// The position weighting method
  /// @brief get weighted combination of position
  ActsVectorD<3> weightPosition() const {
    Vector3D pos = Vector3D(0, 0, 0);
    for (const auto& weightTrack : m_mWeightTracks) {
      pos += weightTrack.first * weightTrack.second.position();
    }
    return pos;
  }

  /// The momentum() weighting method
  /// @brief get weighted combination of momentum
  ActsVectorD<3> weightMomentum() const {
    Vector3D mom = Vector3D(0, 0, 0);
    for (const auto& weightTrack : m_mWeightTracks) {
      mom += weightTrack.first * weightTrack.second.momentum();
    }
    return mom;
  }

  /// The time() weighting method
  /// @brief get weighted combination of time
  double weightTime() const {
    double time = 0.;
    for (const auto& weightTrack : m_mWeightTracks) {
      time += weightTrack.first * weightTrack.second.time();
    }
    return time;
  }

  /// The charge() getter method
  /// @brief get weighted combination of charge
  double weightCharge() const {
    double charge = 0.;
    for (const auto& weightTrack : m_mWeightTracks) {
      charge += weightTrack.first * weightTrack.second.charge();
    }
    return charge;
  }

  /// @note currently get combined full parameter, without covariance
  /// @to do: add combination of covariances
  FullParameterSet weightParameterSet() const {
    std::array<double, 6> pars_array = {{0., 0., 0., 0., 0., 0.}};
    ParVector_t parValues;
    parValues << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4], pars_array[5];
    FullParameterSet parSet(std::nullopt, parValues);
    for (const auto& weightTrack : m_mWeightTracks) {
      parValues += weightTrack.first *
                   weightTrack.second.getParameterSet().getParameters();
    }
    return parSet;
  }

  /// @brief return the size of track map
  size_t size() const { return m_mWeightTracks.size(); }

  /// @brief const getter function of the track map
  const ParameterMapWeightTrack& getTrackList() const {
    return m_mWeightTracks;
  }

  /// @brief non-const getter function of the track map
  ParameterMapWeightTrack& getTrackList() { return m_mWeightTracks; }

  /// @note get combined parameter
  ParVector_t parameters() const { return getParameterSet().getParameters(); }
  /// @brief currently this is in no-use
  /// @to do: add combined covariance
  ParVector_t uncertainty() const { return getParameterSet().getUncertainty(); }

  /// @brief access to the reference surface
  const Surface& referenceSurface() const final { return *m_mSurface; }

  /// @brief access to the the reference frame
  RotationMatrix3D referenceFrame(const GeometryContext& gctx) const final {
    return m_mSurface->transform(gctx).linear();
  }

 private:
  ParameterMapWeightTrack m_mWeightTracks;

  std::shared_ptr<multiSurfaceType> m_mSurface;

  ChargePolicy m_oChargePolicy;    ///< charge policy object distinguishing
                                   /// between charged and neutral tracks
  FullParameterSet m_oParameters;  ///< ParameterSet object holding the
                                   /// parameter values and covariance matrix
  double m_oTime;                  ///< time of the track parametrisation
  ActsVectorD<3> m_vPosition;      ///< 3D vector with global position
  ActsVectorD<3> m_vMomentum;      ///< 3D vector with global momentum
};
}  // namespace Acts
