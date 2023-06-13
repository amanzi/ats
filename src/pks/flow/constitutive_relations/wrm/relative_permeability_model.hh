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
relative permeability as a function of saturation.

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
class RelativePermeabilityModel {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  RelativePermeabilityModel(Teuchos::ParameterList& plist)
    : model_(plist)
  {
    my_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(my_key_);
    s_key_ = Keys::readKey(plist, domain, "saturation", "saturation_liquid");
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    res_ = res[0];
    s_ = deps[1];
  }

  KeyVector getMyKeys() const { return { my_key_ }; }
  KeyVector getDependencies() const { return { s_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    res_(i,0) = model_.k_relative(s_(i,0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    res_(i,0) = model_.d_k_relative(s_(i,0));
  }

 private:
  View_type res_;
  cView_type s_;
  Key my_key_, s_key_;

  WRM_type model_;
};


template <class cView_type, class View_type, class WRM_type>
const std::string RelativePermeabilityModel<cView_type, View_type, WRM_type>::name = "relative permeability " + WRM_type::name;


template<class cView_type, class View_type>
using RelativePermeabilityVanGenuchtenModel = RelativePermeabilityModel<cView_type, View_type, WRMVanGenuchten>;


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
