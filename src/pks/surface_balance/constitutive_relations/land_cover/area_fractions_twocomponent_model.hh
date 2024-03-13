/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the two-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A subgrid model for determining the area fraction of snow vs not snow within a grid cell.
/*!

Uses a simple linear transition to vary between snow-covered and
not-snow-covered ground.

Ordering of the area fractions calculated are: [ground/water, snow].

`"evaluator type`" = `"area fractions, two components`"

.. _area-fractions-twocomponent-evaluator-spec:
.. admonition:: area-fractions-twocomponent-evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
         Mimimum area fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" ``[string]``

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "EvaluatorMultiDOFModelCVByMaterial.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class AreaFractionsTwoComponentModel {
 public:

  static const int n_dependencies = 1;
  static const bool provides_derivatives = false;
  static const int n_dofs = 2;
  static const std::vector<std::string> subfield_names;
  static const std::string eval_type;

  explicit AreaFractionsTwoComponentModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    my_key_ = { Keys::cleanPListName(*plist), tag };
    auto domain = Keys::getDomain(my_key_.first);

    Key domain_snow = Keys::readDomainHint(*plist, domain, "surface", "snow");

    // dependencies
    snow_depth_key_ = Keys::readKeyTag(*plist, domain_snow, "snow depth", "depth", tag);

    min_area_ = plist->get<double>("minimum fractional area [-]", 1.e-5);
    if (min_area_ <= 0.) {
      Errors::Message message(
        "AreaFractionsTwoComponentEvaluator: Minimum fractional area should be > 0.");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    std::string region = model_list.get<std::string>("region");
    land_cover_ = getLandCover(region, model_list, { "snow_transition_depth", });
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];
    sd_ = deps[0];
  }

  void freeViews()
  {
    res_ = View_type();
    sd_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { my_key_ }; }
  KeyTagVector getDependencies() const { return { snow_depth_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    // calculate area of land
    if (sd_(i,0) >= land_cover_.snow_transition_depth) {
      res_(i,1) = 1.;
    } else if (sd_(i,0) <= 0.) {
      res_(i,1) = 0.;
    } else {
      res_(i,1) = sd_(i,0) / land_cover_.snow_transition_depth;
    }

    // if any area is less than eps, give to other
    if (res_(i,1) < min_area_) {
      res_(i,1) = 0.;
    } else if (res_(i,1) > (1 - min_area_)) {
      res_(i,1) = 1.;
    }
    res_(i,0) = 1 - res_(i,1);
  }

  // derivatives not currently provided
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }

 private:
  KeyTag my_key_;
  KeyTag snow_depth_key_;

  View_type res_;
  cView_type sd_;

  LandCover land_cover_;
  double min_area_;
};


template <class cView_type, class View_type>
const std::string AreaFractionsTwoComponentModel<cView_type, View_type>::eval_type =
  "area fractions, two components";

template <class cView_type, class View_type>
const std::vector<std::string> AreaFractionsTwoComponentModel<cView_type, View_type>::subfield_names =
  { "ground_or_water", "snow" };


using AreaFractionsTwoComponentEvaluator = EvaluatorMultiDOFModelCVByMaterial<AreaFractionsTwoComponentModel>;

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
