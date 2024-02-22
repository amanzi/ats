/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Computes the ponded depth given pressure and water density.

.. math::
  h = \frac{p - p_{atm}}{ \rho g }

Note this is only valid for unfrozen water, if variable ice density is used.


type : `"ponded depth`"

.. _ponded-depth-model-spec
.. admonition:: ponded-depth-model-spec

   * `"allow negative ponded depth`" ``[bool]`` **false** If false, sets the
     floor at 0.
                
   KEYS:

   - `"pressure`" **DOMAIN-pressure**
   - `"mass density`" **DOMAIN-mass_density_liquid**

*/

#pragma once

#include "Key.hh"
#include "StateDefs.hh"

namespace Amanzi {

class State;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type>
class PondedDepthModel {
 public:
  static const int n_results = 1;
  static const int n_dependencies = 2;
  static const std::string eval_type; // = "overland pressure water content";

  PondedDepthModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    WC_key_ = Keys::cleanPListName(*plist);
    auto domain = Keys::getDomain(WC_key_);
    pres_key_ = Keys::readKey(*plist, domain, "pressure", "pressure");
    rho_key_ = Keys::readKey(*plist, domain, "mass density", "mass_density_liquid");

    bar_ = plist->get<bool>("allow negative ponded depth", false);
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == n_dependencies);
    AMANZI_ASSERT(res.size() == n_results);
    WC_ = res[0];
    pres_ = deps[0];
    rho_ = deps[1];

    const auto& gravity = s.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    gz_ = -gravity[gravity.dim() - 1];
    p_atm_ = s.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  }

  KeyVector getMyKeys() const { return { WC_key_ }; }
  KeyVector getDependencies() const { return { pres_key_, rho_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int c) const
  {
    if (bar_) {
      WC_(c, 0) = (pres_(c, 0) - p_atm_) / (gz_ * rho_(c, 0));
    } else {
      WC_(c, 0) = pres_(c, 0) < p_atm_ ? 0. : (pres_(c, 0) - p_atm_) / (gz_ * rho_(c, 0));
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int c) const
  {
    if (bar_) {
      WC_(c, 0) = 1.0 / (gz_ * rho_(c, 0));
    } else {
      WC_(c, 0) = pres_(c, 0) < p_atm_ ? 0. : 1.0 / (gz_ * rho_(c, 0));
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int c) const
  {
    if (bar_) {
      WC_(c, 0) = -gz_ * (pres_(c, 0) - p_atm_) / std::pow(gz_ * rho_(c, 0), 2);
    } else {
      WC_(c, 0) =
        pres_(c, 0) < p_atm_ ? 0. : -gz_ * (pres_(c, 0) - p_atm_) / std::pow(gz_ * rho_(c, 0), 2);
    }
  }

 private:
  View_type WC_;
  cView_type pres_, rho_;

  Key WC_key_;
  Key pres_key_, rho_key_;

  bool bar_;
  double gz_;
  double p_atm_;
};

template <class cView_type, class View_type>
const std::string PondedDepthModel<cView_type, View_type>::eval_type = "ponded depth";

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
