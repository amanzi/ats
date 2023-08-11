/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Computes the water content [mols] for a standard Richards equation:

.. math::
  \Theta = n_l s_l \phi |V|

type : `"richards water content`"

.. _richards_water_content_model
.. admonition:: richards-water-content-model-spec

   KEYS:

   - `"density`" **DOMAIN-molar_density_liquid**
   - `"saturation`" **DOMAIN-saturation_liquid**
   - `"porosity`" **DOMAIN-porosity**
   - `"cell volume`" **DOMAIN-cell_volume**

*/

#pragma once

#include "Key.hh"
#include "StateDefs.hh"

namespace Amanzi {

class State;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type>
class RichardsWaterContentModel {
 public:
  static const int n_results = 1;
  static const int n_dependencies = 4;
  static const std::string name; // = "richards water content";

  RichardsWaterContentModel(Teuchos::ParameterList& plist) {
    WC_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(WC_key_);
    nl_key_ = Keys::readKey(plist, domain, "molar density", "molar_density_liquid");
    sl_key_ = Keys::readKey(plist, domain, "saturation", "saturation_liquid");
    phi_key_ = Keys::readKey(plist, domain, "porosity", "porosity");
    cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    AMANZI_ASSERT(deps.size() == n_dependencies);
    AMANZI_ASSERT(res.size() == n_results);
    WC_ = res[0];
    nl_ = deps[0];
    sl_ = deps[1];
    phi_ = deps[2];
    cv_ = deps[3];
  }

  KeyVector getMyKeys() const {
    return { WC_key_ };
  }
  KeyVector getDependencies() const {
    return { nl_key_, sl_key_, phi_key_, cv_key_ };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    WC_(i,0) = nl_(i,0) * sl_(i,0) * phi_(i,0) * cv_(i,0);
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    WC_(i,0) = sl_(i,0) * phi_(i,0) * cv_(i,0);
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    WC_(i,0) = nl_(i,0) * phi_(i,0) * cv_(i,0);
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const {
    WC_(i,0) = nl_(i,0) * sl_(i,0) * cv_(i,0);
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const {
    WC_(i,0) = nl_(i,0) * sl_(i,0) * phi_(i,0);
  }

 private:
  View_type WC_;
  cView_type nl_, sl_, phi_, cv_;

  Key WC_key_;
  Key nl_key_, sl_key_, phi_key_, cv_key_;
};

template <class cView_type, class View_type>
const std::string RichardsWaterContentModel<cView_type, View_type>::name = "richards water content";

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

