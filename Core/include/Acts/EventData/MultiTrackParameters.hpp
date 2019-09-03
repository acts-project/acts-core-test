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
  /// The structure is used in component reduction in GSF
  using ParameterWeightTrack = std::pair<double, ParametersT>;
  using ParameterMapWeightTrack =
      std::multimap<double, ParametersT, std::greater<double>>;

  using multiSurfaceType = typename ParametersT::surfaceType;

  /// @brief standard constructor for bound track parameters of charged particles 
    /// @param[in] weightParList is the list of BoundParameters moved in
	/// @param[in] surface The reference surface the parameters are bound to
  template <typename P = ParametersT, typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters(std::list<std::pair<double, P>>&& weightParList,
                       std::shared_ptr<multiSurfaceType> surface)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_mSurface(std::move(surface)){
		  assert(m_mSurface);
  }

  /// @brief standard constructor for curv track parameters of charged particles
    /// @param[in] weightParList is the list of BoundParameters moved in
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters(
      std::list<
          std::pair<double, SingleCurvilinearTrackParameters<ChargedPolicy>>>&&
          weightParList)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()){
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// @brief standard constructor for bound track parameters of neutral particles
    /// @param[in] weightParList is the list of BoundParameters moved in
	/// @param[in] surface The reference surface the parameters are bound to
  template <typename P = ParametersT, typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters(std::list<std::pair<double, P>>&& weightParList,
                       std::shared_ptr<multiSurfaceType> surface)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()),
        m_mSurface(std::move(surface)){
		  assert(m_mSurface);
  }
  /// @brief standard constructor for curv track parameters of neutral particles
    /// @param[in] weightParList is the list of BoundParameters moved in
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters(
      std::list<
          std::pair<double, SingleCurvilinearTrackParameters<NeutralPolicy>>>&&
          weightParList)
      : TrackParametersBase(),
        m_mWeightTracks(weightParList.begin(), weightParList.end()){
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// @brief copy constructor
  /// @param[in] copy The source parameters
  MultiTrackParameters(
      const MultiTrackParameters<ChargePolicy, ParametersT>& copy) = default;

  /// @brief default move constructor
  /// @param[in] copy The source parameters
  MultiTrackParameters(MultiTrackParameters<ChargePolicy, ParametersT>&& copy) =
      default;

  /// @brief assignment operator for BoundParameters
  MultiTrackParameters<ChargePolicy, ParametersT>& operator=(
      const MultiTrackParameters<ChargePolicy, ParametersT>& rhs) {
    // check for self-assignment
    if (this != &rhs) {
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
      m_mWeightTracks = std::move(rhs.m_mWeightTracks);
      m_mSurface = std::move(rhs.m_mSurface);
    }
    return *this;
  }

  /// @brief default virtual destructor
  ~MultiTrackParameters() override = default;

  /// @brief equality operator
  ///
  bool operator==(const TrackParametersBase& /*unused*/) const final {
    return true;
  }

  /// @brief move a single parameter to the track map for BoundParameters
  /// @param weight the parameter weight
  /// @param para the bound parameter 
  template <typename P = ParametersT>
  void append(double weight, P&& para) {
    // They should be on the same surface
    assert(&para.referenceSurface() == m_mSurface.get());
    m_mWeightTracks.insert(std::make_pair(weight, std::move(para)));
  }

  /// @brief move a single parameter to the track map for CurvilinearParameters
  /// @param weight the parameter weight
  /// @param para the curvilinear parameter 
  void append(double weight,
              SingleCurvilinearTrackParameters<ChargePolicy>&& para) {
    m_mWeightTracks.insert(std::make_pair(weight, std::move(para)));
    m_mSurface = Surface::makeShared<PlaneSurface>(position(), momentum());
  }

  /// getter functions for position 
  ActsVectorD<3> position() const final {
    Vector3D pos = Vector3D(0, 0, 0);
    for (const auto& weightTrack : m_mWeightTracks) {
      pos += weightTrack.first * weightTrack.second.position();
    }
    return pos;
  }
  /// getter functions for momentum
  ActsVectorD<3> momentum() const final { 
    Vector3D mom = Vector3D(0, 0, 0);
    for (const auto& weightTrack : m_mWeightTracks) {
      mom += weightTrack.first * weightTrack.second.momentum();
    }
    return mom;
  }
  /// getter functions for time
  double time() const final {
    double time = 0.;
    for (const auto& weightTrack : m_mWeightTracks) {
      time += weightTrack.first * weightTrack.second.time();
    }
    return time;
  }
  /// getter functions for charge
  double charge() const final { 
	return  m_mWeightTracks.begin()->second.charge();
  }

  /// getter functions for parset
  const FullParameterSet& getParameterSet() const final {
	// Combine parameters 
    std::array<double, 6> pars_array = {{0., 0., 0., 0., 0., 0.}};
    ParVector_t parsCombine;
    parsCombine<< pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4], pars_array[5];
    for (const auto& weightTrack : m_mWeightTracks) {
      parsCombine += weightTrack.first *
                   weightTrack.second.getParameterSet().getParameters();
    }

	// Combine covariances
	CovMatrix_t covPart_1 = CovMatrix_t::Zero();
	CovMatrix_t covPart_2 = CovMatrix_t::Zero();
	typename ParameterMapWeightTrack::const_iterator iter = m_mWeightTracks.begin();
	double weightSum = 0.;
	for ( ;iter != m_mWeightTracks.end(); iter++){
			weightSum += (*iter).first;
			covPart_1 += (*iter).first * (*(*iter).second.covariance());
	typename ParameterMapWeightTrack::const_iterator iterRemaining = iter;
	for ( ;iterRemaining != m_mWeightTracks.end(); iterRemaining++ ){
	  if( iterRemaining == iter ) continue;
	  auto paraDiff = (*iter).second.parameters() - (*iterRemaining).second.parameters();
	  auto unity  = paraDiff * paraDiff.transpose();
	  covPart_2 = (*iter).first * (*iterRemaining).first * unity;
	}
	}
	CovMatrix_t covCombine = covPart_1/weightSum + covPart_2/(weightSum * weightSum);

    FullParameterSet parSet(std::move(covCombine), parsCombine);
		
	// return to a const& is wrong
	// @to be improved
    return std::move(parSet);
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
  ParVector_t parameters() const { 
	return getParameterSet().getParameters(); 
  }
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

};
}  // namespace Acts
