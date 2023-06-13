/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Uses WRM class instances to compute WRMs.
/*!

This is a lightweight wrapper around things that can behave as WRMs, evaluating
saturation as a function of capillary pressure.

*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"

#include "wrm_van_genuchten.hh"

namespace Amanzi {

class State;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type, class WRM_type>
class WRMModel {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  WRMModel(Teuchos::ParameterList& plist)
    : model_(plist)
  {
    my_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(my_key_);
    pc_key_ = Keys::readKey(plist, domain, "capillary pressure", "capillary_pressure_liq_atm");
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    res_ = res[0];
    pc_ = deps[1];
  }

  KeyVector getMyKeys() const { return { my_key_ }; }
  KeyVector getDependencies() const { return { pc_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    res_(i,0) = model_.saturation(pc_(i,0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    res_(i,0) = model_.d_saturation(pc_(i,0));
  }

 private:
  View_type res_;
  cView_type pc_;
  Key my_key_, pc_key_;

  WRM_type model_;
};


template <class cView_type, class View_type, class WRM_type>
const std::string WRMModel<cView_type, View_type, WRM_type>::name = "wrm " + WRM_type::name;


template<class cView_type, class View_type>
using WRMVanGenuchtenModel = WRMModel<cView_type, View_type, WRMVanGenuchten>;


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
