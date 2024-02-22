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
  static const std::string eval_type;

  WRMModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : model_(plist->sublist("model parameters"))
  {
    // my keys are for saturation, note that order matters, liquid -> gas
    Key akey = Keys::cleanPListName(*plist);
    Key domain_name = Keys::getDomain(akey);

    std::size_t liq_pos = akey.find("liquid");
    std::size_t gas_pos = akey.find("gas");
    if (liq_pos != std::string::npos) {
      sl_key_ = akey;
      Key otherkey = akey.substr(0, liq_pos) + "gas" + akey.substr(liq_pos + 6);
      sg_key_ = Keys::readKey(*plist, domain_name, "other saturation", otherkey);

    } else if (gas_pos != std::string::npos) {
      sl_key_ = akey.substr(0, gas_pos) + "liquid" + akey.substr(gas_pos + 3);
      sl_key_ = Keys::readKey(*plist, domain_name, "saturation", sl_key_);
      sg_key_ = akey;

    } else {
      sl_key_ = Keys::readKey(*plist, domain_name, "saturation");
      sg_key_ = Keys::readKey(*plist, domain_name, "other saturation");
    }

    pc_key_ =
      Keys::readKey(*plist, domain_name, "capillary pressure", "capillary_pressure_gas_liq");
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    sl_ = res[0];
    sg_ = res[1];
    pc_ = deps[0];
  }

  KeyVector getMyKeys() const { return { sl_key_, sg_key_ }; }
  KeyVector getDependencies() const { return { pc_key_ }; }
  WRM_type& getModel() { return model_; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    sl_(i, 0) = model_.saturation(pc_(i, 0));
    sg_(i, 0) = 1 - sl_(i, 0);
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    sl_(i, 0) = model_.d_saturation(pc_(i, 0));
    sg_(i, 0) = -sl_(i, 0);
  }

 private:
  View_type sl_, sg_;
  cView_type pc_;
  Key pc_key_;
  Key sl_key_, sg_key_;

  WRM_type model_;
};


template <class cView_type, class View_type, class WRM_type>
const std::string WRMModel<cView_type, View_type, WRM_type>::eval_type = "wrm " + WRM_type::eval_type;


template <class cView_type, class View_type>
using WRMVanGenuchtenModel = WRMModel<cView_type, View_type, WRMVanGenuchten>;


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
