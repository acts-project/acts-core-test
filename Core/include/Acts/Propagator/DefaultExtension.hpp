// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// @brief Default evaluater of the k_i's and elements of the transport matrix
/// D of the RKN4 stepping. This is a pure implementation by textbook.
struct DefaultExtension
{
  /// @brief Default constructor
  DefaultExtension() = default;

  /// Local store for q/p in SI units
  double qop = 0.;

  /// Local store for conversion of momentum from SI to natural units
  const double m_conv = units::SI2Nat<units::MOMENTUM>(1);

  /// @brief Control function if the step evaluation would be valid
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  /// @return Boolean flag if the step would be valid
  template <typename stepper_state_t>
  bool
  validExtensionForStep(const stepper_state_t& state) const
  {
    // Check existence of a volume with material
    if (state.volume && (*state.volume) && (*state.volume)->material())
      return false;
    return true;
  }

  /// @brief Evaluater of the k_i's of the RKN4. For the case of i = 0 this
  /// step sets up qop, too.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  /// @param [out] knew Next k_i that is evaluated
  /// @param [in] bField B-Field at the evaluation position
  /// @param [in] i Index of the k_i
  /// @param [in] h Step size (= 0. ^ 0.5 * StepSize ^ StepSize)
  /// @param [in] kprev Evaluated k_{i - 1}
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t>
  bool
  k(const stepper_state_t& state,
    Vector3D&              knew,
    const Vector3D&        bField,
    const int              i     = 0,
    const double           h     = 0.,
    const Vector3D&        kprev = Vector3D())
  {
    // First step does not rely on previous data
    if (i == 0) {
      // Store qop, it is always used if valid
      qop = state.q / units::Nat2SI<units::MOMENTUM>(state.p);

      // Evaluate the k_i
      knew = qop * state.dir.cross(bField);
    } else {
      knew = qop * (state.dir + h * kprev).cross(bField);
    }
    return true;
  }

  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Since the textbook does not deliver further vetos,
  /// this is a dummy function.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t>
  bool
  finalize(stepper_state_t& /*unused*/, const double /*unused*/) const
  {
    return true;
  }

  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Since the textbook does not deliver further vetos,
  /// this is just for the evaluation of the transport matrix.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @tparam stepper_data_t Type of the data collected in the step
  /// @param [in] state State of the stepper
  /// @param [in] h Step size
  /// @param [in] data Data of B-field and k_i's
  /// @param [out] D Transport matrix
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t, typename stepper_data_t>
  bool
  finalize(stepper_state_t&      state,
           const double          h,
           const stepper_data_t& data,
           ActsMatrixD<7, 7>& D) const
  {
    return transportMatrix(state.dir, h, data, D);
  }

protected:
  /// @brief Evaluates the transport matrix D for the jacobian
  ///
  /// @param [in] dir Direction of the particle
  /// @param [in] h Step size
  /// @param [in] sd Data of B-field and k_i's
  /// @param [out] D Transport matrix
  /// @return Boolean flag if evaluation is valid
  template <typename stepper_data_t>
  bool
  transportMatrix(const Vector3D&       dir,
                  const double          h,
                  const stepper_data_t& sd,
                  ActsMatrixD<7, 7>& D) const
  {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by dFdT and dGdT. The second is given by the top 3
    /// lines of the rightmost column. This is calculated by dFdL and dGdL.
    /// The remaining non-zero term is calculated directly. The naming of the
    /// variables is explained in eq. 11 and are directly related to the
    /// initial problem in eq. 7.
    /// The evaluation is based by propagating the parameters T and lambda as
    /// given in eq. 16 and evaluating the derivations for matrix A.

    D = ActsMatrixD<7, 7>::Identity();

    double half_h = h * 0.5;
    // This sets the reference to the sub matrices
    // dFdx is already initialised as (3x3) zero
    auto dFdT = D.block<3, 3>(0, 3);
    auto dFdL = D.block<3, 1>(0, 6);
    // dGdx is already initialised as (3x3) identity
    auto dGdT = D.block<3, 3>(3, 3);
    auto dGdL = D.block<3, 1>(3, 6);

    ActsMatrixD<3, 3> dk1dT = ActsMatrixD<3, 3>::Zero();
    ActsMatrixD<3, 3> dk2dT = ActsMatrixD<3, 3>::Identity();
    ActsMatrixD<3, 3> dk3dT = ActsMatrixD<3, 3>::Identity();
    ActsMatrixD<3, 3> dk4dT = ActsMatrixD<3, 3>::Identity();

    ActsVectorD<3> dk1dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk2dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk3dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk4dL = ActsVectorD<3>::Zero();

    // For the case without energy loss
    dk1dL = dir.cross(sd.B_first);
    dk2dL = (dir + half_h * sd.k1).cross(sd.B_middle)
        + qop * half_h * dk1dL.cross(sd.B_middle);
    dk3dL = (dir + half_h * sd.k2).cross(sd.B_middle)
        + qop * half_h * dk2dL.cross(sd.B_middle);
    dk4dL
        = (dir + h * sd.k3).cross(sd.B_last) + qop * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= qop;

    dk2dT += half_h * dk1dT;
    dk2dT *= cross(dk2dT, sd.B_middle);
    dk2dT *= qop;

    dk3dT += half_h * dk2dT;
    dk3dT *= cross(dk3dT, sd.B_middle);
    dk3dT *= qop;

    dk4dT += h * dk3dT;
    dk4dT *= cross(dk4dT, sd.B_last);
    dk4dT *= qop;

    dFdT.setIdentity();
    dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = m_conv * h * h / 6 * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

    dGdL = m_conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

    return true;
  }

  /// @brief Calculates column-wise cross products of a matrix and a vector and
  /// stores the result column-wise in a matrix.
  ///
  /// @param [in] m Matrix that will be used for cross products
  /// @param [in] v Vector for cross products
  /// @return Constructed matrix
  ActsMatrixD<3, 3>
  cross(const ActsMatrixD<3, 3>& m, const Vector3D& v) const
  {
    ActsMatrixD<3, 3> r;
    r.col(0) = m.col(0).cross(v);
    r.col(1) = m.col(1).cross(v);
    r.col(2) = m.col(2).cross(v);

    return r;
  }
};
}  // namespace Acts