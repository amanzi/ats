/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Computes the water content [mols] as a function of surface pressure.

.. math::
  \Theta = \frac{p - p_{atm}}{ M g } |V|

type : `"overland pressure water content`"

.. _overland-pressure-water-content-model-spec
.. admonition:: overland-pressure-water-content-model-spec

   * `"molar mass [kg]`" ``[double]`` **0.0180153** mass of one mol of fluid
   * `"allow negative water content`" ``[bool]`` **false** If false, sets the
     floor at 0.
   * `"water content rollover [Pa]`" ``[double]`` **0** If non-zero, water
     content is splined to have zero derivative (with respect to p) at 0, then
     grow to have derivative 1 at p = this value.
                
   KEYS:

   - `"pressure`" **DOMAIN-pressure**
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
class OverlandPressureWaterContentModel {
 public:
  static const int n_results = 1;
  static const int n_dependencies = 2;
  static const std::string eval_type; // = "overland pressure water content";

  OverlandPressureWaterContentModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    WC_key_ = Keys::cleanPListName(*plist);
    auto domain = Keys::getDomain(WC_key_);
    pres_key_ = Keys::readKey(*plist, domain, "pressure", "pressure");
    cv_key_ = Keys::readKey(*plist, domain, "cell volume", "cell_volume");

    M_ = plist->get<double>("molar mass", 0.0180153);
    bar_ = plist->get<bool>("allow negative water content", false);
    rollover_ = plist->get<double>("water content rollover", 0.);
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == n_dependencies);
    AMANZI_ASSERT(res.size() == n_results);
    WC_ = res[0];
    pres_ = deps[0];
    cv_ = deps[1];

    const auto& gravity = s.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    Mgz_ = -gravity[gravity.dim() - 1] * M_;
    p_atm_ = s.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  }

  KeyVector getMyKeys() const { return { WC_key_ }; }
  KeyVector getDependencies() const { return { pres_key_, cv_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int c) const
  {
    if (bar_) {
      WC_(c, 0) = cv_(c, 0) * (pres_(c, 0) - p_atm_) / Mgz_;
    } else if (rollover_ > 0.) {
      double dp = pres_(c, 0) - p_atm_;
      double dp_eff = dp < 0.        ? 0. :
                      dp < rollover_ ? dp * dp / (2 * rollover_) :
                                       dp - rollover_ / 2.;
      WC_(c, 0) = cv_(c, 0) * dp_eff / Mgz_;
    } else {
      WC_(c, 0) = pres_(c, 0) < p_atm_ ? 0. : cv_(c, 0) * (pres_(c, 0) - p_atm_) / Mgz_;
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int c) const
  {
    if (bar_) {
      WC_(c, 0) = cv_(c, 0) / Mgz_;
    } else if (rollover_ > 0.) {
      double dp = pres_(c, 0) - p_atm_;
      double ddp_eff = dp < 0. ? 0. : dp < rollover_ ? dp / rollover_ : 1.;
      WC_(c, 0) = cv_(c, 0) * ddp_eff / Mgz_;
    } else {
      WC_(c, 0) = pres_(c, 0) < p_atm_ ? 0. : cv_(c, 0) / Mgz_;
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int c) const
  {
    if (bar_) {
      WC_(c, 0) = (pres_(c, 0) - p_atm_) / Mgz_;
    } else if (rollover_ > 0.) {
      double dp = pres_(c, 0) - p_atm_;
      double dp_eff = dp < 0.        ? 0. :
                      dp < rollover_ ? dp * dp / (2 * rollover_) :
                                       dp - rollover_ / 2.;
      WC_(c, 0) = dp_eff / Mgz_;
    } else {
      WC_(c, 0) = pres_(c, 0) < p_atm_ ? 0. : (pres_(c, 0) - p_atm_) / Mgz_;
    }
  }

 private:
  View_type WC_;
  cView_type pres_, cv_;

  Key WC_key_;
  Key pres_key_, cv_key_;

  double M_;
  bool bar_;
  double rollover_;

  double Mgz_;
  double p_atm_;
};

template <class cView_type, class View_type>
const std::string OverlandPressureWaterContentModel<cView_type, View_type>::eval_type =
  "overland pressure water content";

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
