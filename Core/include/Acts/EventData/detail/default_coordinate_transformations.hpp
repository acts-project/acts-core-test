// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

static ParVector_t
global2curvilinear(const ActsVectorD<3>& /*pos*/,
                   const ActsVectorD<3>& mom,
                   double                charge)
{
  using VectorHelpers::phi;
  using VectorHelpers::theta;
  ParVector_t parameters = ParVector_t::Zero();
  parameters(2)          = phi(mom);
  parameters(3)          = theta(mom);
  parameters(4) = ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm();

  return parameters;
}

static ParVector_t
global2parameters(const ActsVectorD<3>& pos,
                  const ActsVectorD<3>& mom,
                  double                charge,
                  const Surface&        s)
{
  using VectorHelpers::phi;
  using VectorHelpers::theta;
  ActsVectorD<2> localPosition;
  s.globalToLocal(pos, mom, localPosition);
  ParVector_t result = ParVector_t::Zero();
  result(0)          = localPosition(0);
  result(1)          = localPosition(1);
  result(2)          = phi(mom);
  result(3)          = theta(mom);
  result(4)          = ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm();

  return result;
}