/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*

Capillary pressure as a function of absolute pressure, given by:

.. math::
  pc = p_atm - p

type : `"capillary pressure, atmospheric gas over liquid`"

.. _capillary-pressure-liquid-atm-model-spec
.. admonition:: capillary-pressure-liquid-atm-model-spec

   KEYS:
   - `"pressure`"

*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"

namespace Amanzi {

class State;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type>
class CapillaryPressureLiquidAtmModel {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  CapillaryPressureLiquidAtmModel(Teuchos::ParameterList& plist) {
    my_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(my_key_);
    p_key_ = Keys::readKey(plist, domain, "pressure", "pressure");
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    p_atm_ = s.Get<double>("atmospheric_pressure", Tags::DEFAULT);
    res_ = res[0];
    p_ = deps[0];
  }

  KeyVector getMyKeys() const { return { my_key_ }; }
  KeyVector getDependencies() const { return { p_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    res_(i,0) = p_atm_ - p_(i,0);
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    res_(i,0) = -1.;
  }

 private:
  View_type res_;
  cView_type p_;
  double p_atm_;

  Key my_key_;
  Key p_key_;
};

template <class cView_type, class View_type>
const std::string CapillaryPressureLiquidAtmModel<cView_type, View_type>::name = "capillary pressure, atmospheric gas over liquid";

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

