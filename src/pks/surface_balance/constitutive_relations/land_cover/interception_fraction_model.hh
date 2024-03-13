/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Partitions water between interception and throughfall, combining throughfall with drainage.
/*!

Based on CLM 4.5 and Lawrence et al 2007, interception is given by:

.. math::
   I = (P_{rain} + P_{snow}) * \alpha * (1 - exp(-.5(LAI)))

Throughfall is given by:

.. math::
   T = (P_{rain} + P_{snow}) - I

Drainage is provided as input here, as a total drainage from the canopy.  The
phase of this drainage is assumed to match the phase of the precipitation.  So
if it is raining, drainage is rain, while if it is 50/50 rain and snow,
drainage is also 50/50 rain and snow.  If total precipitation is 0, then
drainage is partitioned by air temperature (above 0C --> all rain, otherwise
all snow).  This evaluator partitions the drainage and sums it with throughfall
to compute the total source, in each phase, to the layer below the canopy (snow
and/or ground surface).

.. _interception-fraction-model-spec:
.. admonition:: interception-fraction-model-spec

   * `"region`" ``[string]`` Region on which this is applied.
   * `"leaf area interception fraction [-]`" ``[double]`` **0.25** The alpha
      term, this describes the fraction of leaf area that intercepts water.

   MY KEYS:
   - `"interception`" **DOMAIN-interception**
   - `"throughfall and drainage rain`" **DOMAIN-throughfall_drainage_rain**
   - `"throughfall and drainage snow`" **DOMAIN-throughfall_drainage_snow**

   KEYS:
   - `"area index`" **DOMAIN-area_index**
   - `"precipitation rain`" **DOMAIN_SURFACE-precipitation_rain**
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation**
   - `"drainage`" **DOMAIN-drainage**
   - `"air temperature`" **DOMAIN_SURFACE-air_temperature**

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
class InterceptionFractionModel {
 public:

  static const int n_dependencies = 5;
  static const bool provides_derivatives = true;
  static const std::string eval_type;

  explicit InterceptionFractionModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Key akey = Keys::cleanPListName(*plist);
    Key domain = Keys::getDomain(akey);
    Tag tag(plist->get<std::string>("tag"));

    akey = Keys::getVarName(akey);
    Key domain_surf = Keys::readDomainHint(*plist, domain, "canopy", "surface");
    Key domain_snow = Keys::readDomainHint(*plist, domain, "canopy", "snow");

    // my keys
    Key interception_key = Keys::in(akey, "interception") ? akey : "interception";
    interception_key_ = Keys::readKeyTag(*plist, domain, "interception", interception_key, tag);

    Key throughfall_rain_key =
      (Keys::in(akey, "rain") && Keys::in(akey, "throughfall")) ? akey : "throughfall_drainage_rain";
    throughfall_rain_key_ =
      Keys::readKeyTag(*plist, domain, "throughfall and drainage rain", throughfall_rain_key, tag);

    Key throughfall_snow_key =
      (Keys::in(akey, "snow") && Keys::in(akey, "throughfall")) ? akey : "throughfall_drainage_snow";
    throughfall_snow_key_ =
      Keys::readKeyTag(*plist, domain, "throughfall and drainage snow", throughfall_snow_key, tag);

    // - pull Keys from plist
    // dependency: surface-area_index
    ai_key_ = Keys::readKeyTag(*plist, domain, "area index", "area_index", tag);
    rain_key_ = Keys::readKeyTag(*plist, domain_surf, "precipitation rain", "precipitation_rain", tag);
    snow_key_ = Keys::readKeyTag(*plist, domain_snow, "precipitation snow", "precipitation", tag);
    drainage_key_ = Keys::readKeyTag(*plist, domain, "drainage", "drainage", tag);
    air_temp_key_ = Keys::readKeyTag(*plist, domain_surf, "air temperature", "air_temperature", tag);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    alpha_ = model_list.get<double>("leaf area interception fraction [-]", 0.25);
    if (alpha_ < 0 || alpha_ > 1) {
      Errors::Message msg;
      msg << "InterceptionFraction: invalid \"leaf area interception fraction [-]\", must be in "
        "[0,1] (provided: "
          << alpha_ << ")";
      Exceptions::amanzi_throw(msg);
    }
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    interception_ = res[0];
    throughfall_rain_ = res[1];
    throughfall_snow_ = res[2];

    ai_ = deps[0];
    rain_ = deps[1];
    snow_ = deps[2];
    drainage_ = deps[3];
    air_temp_ = deps[4];
  }

  void freeViews()
  {
    interception_ = View_type();
    throughfall_rain_ = View_type();
    throughfall_snow_ = View_type();

    ai_ = cView_type();
    rain_ = cView_type();
    snow_ = cView_type();
    drainage_ = cView_type();
    air_temp_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { interception_key_, throughfall_rain_key_, throughfall_snow_key_ }; }
  KeyTagVector getDependencies() const { return { ai_key_, rain_key_, snow_key_, drainage_key_, air_temp_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double coef = alpha_ * (1 - exp(-0.5 * ai_(i,0)));
    double total_precip = rain_(i,0) + snow_(i,0);
    interception_(i,0) = total_precip * coef;

    double frac_r =
      total_precip > 0 ? rain_(i,0) / total_precip : (air_temp_(i,0) > 273.15 ? 1 : 0);
    throughfall_rain_(i,0) = (1 - coef) * rain_(i,0) + frac_r * drainage_(i,0);
    throughfall_snow_(i,0) = (1 - coef) * snow_(i,0) + (1 - frac_r) * drainage_(i,0);
  }

  // d/d_ai
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    assert(false);
  }

  // d/d_rain
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    assert(false);
  }

  // d/d_snow
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const {
    assert(false);
  }

  // d/d_drainage
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const
  {
    interception_(i,0) = 0.;

    double total_precip = rain_(i,0) + snow_(i,0);
    double frac_r =
      total_precip > 0 ? rain_(i,0) / total_precip : (air_temp_(i,0) > 273.15 ? 1 : 0);
    throughfall_rain_(i,0) = frac_r;
    throughfall_snow_(i,0) = 1 - frac_r;
  }

  // d/d_air_temp
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<4>, const int i) const {
    assert(false);
  }

 private:
  View_type interception_, throughfall_rain_, throughfall_snow_;
  cView_type ai_, rain_, snow_, drainage_, air_temp_;

  KeyTag interception_key_, throughfall_snow_key_, throughfall_rain_key_;
  KeyTag ai_key_, rain_key_, snow_key_, drainage_key_, air_temp_key_;

  double alpha_;
};


template <class cView_type, class View_type>
const std::string InterceptionFractionModel<cView_type, View_type>::eval_type =
  "interception fraction";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
