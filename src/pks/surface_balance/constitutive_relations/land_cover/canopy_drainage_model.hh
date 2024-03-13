/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Drainage rate from the canopy to the lower layers.
/*!

A simple model based on relaxation from current water content to a saturated water content.

.. code::

          |
          | source
          V
         /   \
      I /     \
       V       |
   --Theta--    | T
       ^       |
       | D     |
       V       V
   -- -- -- -- -- -- --


This is the model for drainage D.

Drainage is given by:

.. math::
   D = max(0, \frac{(\Theta - \Theta_sat)}{\tau})

.. _drainage-evaluator-spec:
.. admonition:: drainage-evaluator-spec

   * `"drainage timescale [s]`" ``[double]`` **864** Timescale over which drainage occurs.
   * `"saturated specific water content [m^3 H2O / m^2 leaf area]`" ``[double]`` **1e-4**
      The thickness of the wetting surface -- determines maximum canopy water storage.\

   KEYS:
   - "area index"
   - "water equivalent"

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
class CanopyDrainageModel {
 public:

  static const int n_dependencies = 2;
  static const bool provides_derivatives = true;
  static const std::string eval_type;

  explicit CanopyDrainageModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Key akey = Keys::cleanPListName(*plist);
    Key domain = Keys::getDomain(akey);
    Tag tag(plist->get<std::string>("tag"));

    akey = Keys::getVarName(akey);

    // my keys
    Key drainage_key = Keys::in(akey, "drainage") ? akey : "drainage";
    drainage_key_ = Keys::readKeyTag(*plist, domain, "drainage", drainage_key, tag);

    Key fracwet_key = Keys::in(akey, "fracwet") ? akey : "fracwet";
    fracwet_key_ = Keys::readKeyTag(*plist, domain, "fraction wet", fracwet_key, tag);

    // dependencies
    ai_key_ = Keys::readKeyTag(*plist, domain, "area index", "area_index", tag);
    wc_key_ = Keys::readKeyTag(*plist, domain, "water equivalent", "water_equivalent", tag);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    tau_ = model_list.get<double>("drainage timescale [s]", 864);
    if (tau_ <= 0) {
      Errors::Message msg("CanopyDrainageModel: invalid \"drainage timescale [s]\", must be positive.");
      Exceptions::amanzi_throw(msg);
    }

    wc_sat_ = model_list.get<double>("saturated specific water content [m^3 H2O / m^2 leaf area]", 1.e-4);
    if (wc_sat_ < 0) {
      Errors::Message msg("\"saturated specific water content\" must be greater than 0.");
      Exceptions::amanzi_throw(msg);
    }
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    drainage_ = res[0];
    fracwet_ = res[1];

    ai_ = deps[0];
    wc_ = deps[1];
  }

  void freeViews()
  {
    drainage_ = View_type();
    fracwet_ = View_type();
    ai_ = cView_type();
    wc_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { drainage_key_, fracwet_key_ }; }
  KeyTagVector getDependencies() const { return { ai_key_, wc_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double wc_cell = fmax(wc_(i,0), 0.);
    double ai_cell = fmax(ai_(i,0), 0.);

    // convert from m^3 H20/ m^2 leaf area to m^3 H20 / m^2 cell area
    double wc_cell_sat = wc_sat_ * ai_cell;
    fracwet_(i,0) = wc_cell_sat > 0. ? wc_cell / wc_cell_sat : 0.;

    // must be in [0,1] -- note that wc_cell can be > wc_cell_sat
    fracwet_(i,0) = fmax(fmin(fracwet_(i,0), 1.0), 0.0);

    if (wc_cell > wc_cell_sat) {
      //  is oversaturated and draining
      // NOTE: should this actually be:
      // res_drainage_c(i,0) = (wc_cell - wc_cell_sat) / ai(i,0) / tau_;
      // to make it proportional in units of m^3 H20 per m^2 leaf area
      // but then we would have to multiply by ai to use it as a source
      drainage_(i,0) = (wc_cell - wc_cell_sat) / tau_;
    } else {
      drainage_(i,0) = 0.;
    }
  }

  // d/d_ai
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    assert(false);
  }

  // d/d_wc
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    double wc_cell = fmax(wc_(i,0), 0.);
    double wc_cell_sat = wc_sat_ * ai_(i,0);
    fracwet_(i,0) = wc_cell_sat > 0. ? 1 / wc_cell_sat : 0;
    if (wc_cell > wc_cell_sat) {
      //  is oversaturated and draining
      drainage_(i,0) = 1.0 / tau_;
    } else {
      drainage_(i,0) = 0.;
    }
  }

 private:
  View_type drainage_, fracwet_;
  cView_type ai_, wc_;

  KeyTag drainage_key_, fracwet_key_;
  KeyTag ai_key_, wc_key_;

  double tau_, wc_sat_;
};


template <class cView_type, class View_type>
const std::string CanopyDrainageModel<cView_type, View_type>::eval_type =
  "canopy drainage";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
