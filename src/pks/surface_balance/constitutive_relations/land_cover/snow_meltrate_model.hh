/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates snow melt via USDA - Natural Resources Conservation Service model
/*!

From:  National Engineering Handbook (NEH) part 630, Chapter 11

Snow melt rate is given by:

.. math::
   SM = H(T_{snow}^{expected} - 273.15) R

where :math:`R` is the snow melt rate per degree-day and
:math:`T_{snow}^{expected}` is the expected snow temperature, which is
typically given by :math:`T_{snow}^{expected} = T_{air} - \Delta`, where
:math:`\Delta` is the expected air-snow temperature difference [C] (note this
is NOT prescribed here -- the user must supply the expected snow temperature
via an evalutor).

Note that the Heaviside function is used to ensure this is only active when the
expected snow temperature is above 0 C.

Then, a linear transition factor is applied to ensure that this goes to zero as
the snow SWE goes to zero.  That factor is 1 at the snow transition depth, and
0 when snow SWE is 0.  This uses LandCover for the snow_ground_transition
parameter.

.. _snow-meltrate-evaluator-spec:
.. admonition:: snow-meltrate-evaluator-spec

   * `"snow melt rate [mm day^-1 C^-1]`" ``[double]`` **2.74**
     the melt rate per degree-day, above 0 C, e.g. :math:`R` above.

   KEYS:

   - `"snow water equivalent`" **DOMAIN-water_equivalent**
   - `"expected snow temperature`"  **DOMAIN-expected_temperature**

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class SnowMeltRateModel {
 public:
  static const int n_dependencies = 2;
  static const bool provides_derivatives = true;
  static const std::string eval_type;

  explicit SnowMeltRateModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    my_key_ = { Keys::cleanPListName(*plist), tag };
    auto domain = Keys::getDomain(my_key_.first);

    exp_temp_key_ =
      Keys::readKeyTag(*plist, domain, "expected snow temperature", "expected_temperature", tag);
    swe_key_ = Keys::readKeyTag(*plist, domain, "snow water equivalent", "water_equivalent", tag);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    melt_rate_ = model_list.get<double>("snow melt rate [mm day^-1 C^-1]");
    melt_rate_ *= 0.001 / 86400.; // convert mm/day to m/s

    if (melt_rate_ < 0) {
      Errors::Message msg;
      msg << "SnowMeltRateModel: invalid \"snow melt rate [mm day^-1 C^-1]\", must be positive.";
      Exceptions::amanzi_throw(msg);
    }

    std::string region = model_list.get<std::string>("region");
    land_cover_ = getLandCover(region, model_list, { "snow_transition_depth" });
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];
    exp_temp_ = deps[0];
    swe_ = deps[1];
  }

  void freeViews()
  {
    res_ = View_type();
    exp_temp_ = cView_type();
    swe_ = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return {
      my_key_,
    };
  }
  KeyTagVector getDependencies() const { return { exp_temp_key_, swe_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    if (exp_temp_(i, 0) > 273.15) {
      res_(i, 0) = melt_rate_ * (exp_temp_(i, 0) - 273.15);

      if (swe_(i, 0) < land_cover_.snow_transition_depth) {
        res_(i, 0) *= Kokkos::max(0., swe_(i, 0) / land_cover_.snow_transition_depth);
      }

    } else {
      res_(i, 0) = 0.0;
    }
  }

  // d/d_ai
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    if (exp_temp_(i, 0) > 273.15) {
      res_(i, 0) = melt_rate_;
      if (swe_(i, 0) < land_cover_.snow_transition_depth) {
        res_(i, 0) *= Kokkos::max(0., swe_(i, 0) / land_cover_.snow_transition_depth);
      }
    } else {
      res_(i, 0) = 0.0;
    }
  }

  // d/d_rain
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    if (swe_(i, 0) < land_cover_.snow_transition_depth && exp_temp_(i, 0) > 273.15) {
      res_(i, 0) = melt_rate_ * (exp_temp_(i, 0) - 273.15) / land_cover_.snow_transition_depth;
    } else {
      res_(i, 0) = 0.0;
    }
  }

 private:
  View_type res_;
  cView_type exp_temp_, swe_;

  KeyTag my_key_;
  KeyTag exp_temp_key_, swe_key_;

  LandCover land_cover_;
  double melt_rate_;
};


template <class cView_type, class View_type>
const std::string SnowMeltRateModel<cView_type, View_type>::eval_type = "snow melt rate";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
