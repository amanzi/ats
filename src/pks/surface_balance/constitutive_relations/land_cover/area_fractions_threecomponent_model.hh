/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell.
/*!

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground, water, snow].

`"evaluator type`" = `"area fractions, three components`"

.. _area-fractions-threecomponent-evaluator-spec:
.. admonition:: area-fractions-threecomponent-evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
      Mimimum area fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" **DOMAIN_SNOW-depth**
   - `"ponded depth`" **DOMAIN-ponded_depth**

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.
   - `"water transition height [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the water gradient.

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
class AreaFractionsThreeComponentModel {
 public:

  static const int n_dependencies = 2;
  static const bool provides_derivatives = false;
  static const int n_dofs = 3;
  static const std::vector<std::string> subfield_names;
  static const std::string eval_type;

  explicit AreaFractionsThreeComponentModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    my_key_ = { Keys::cleanPListName(*plist), tag };
    auto domain = Keys::getDomain(my_key_.first);

    Key domain_snow = Keys::readDomainHint(*plist, domain, "surface", "snow");

    // dependencies
    ponded_depth_key_ = Keys::readKeyTag(*plist, domain, "ponded depth", "ponded_depth", tag);
    snow_depth_key_ = Keys::readKeyTag(*plist, domain_snow, "snow depth", "depth", tag);

    min_area_ = plist->get<double>("minimum fractional area [-]", 1.e-5);
    if (min_area_ <= 0.) {
      Errors::Message message(
        "AreaFractionsThreeComponentEvaluator: Minimum fractional area should be > 0.");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    std::string region = model_list.get<std::string>("region");
    land_cover_ = getLandCover(region, model_list, { "snow_transition_depth", "water_transition_depth" });
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];
    pd_ = deps[0];
    sd_ = deps[1];
  }

  void freeViews()
  {
    res_ = View_type();
    pd_ = cView_type();
    sd_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { my_key_ }; }
  KeyTagVector getDependencies() const { return { ponded_depth_key_, snow_depth_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    // calculate area of land
    if (sd_(i,0) >= land_cover_.snow_transition_depth) {
      res_(i,1) = 0.;
      res_(i,2) = 1.;
    } else {
      if (sd_(i,0) <= 0.) {
        res_(i,2) = 0.;
      } else {
        res_(i,2) = sd_(i,0) / land_cover_.snow_transition_depth;
      }

      // snow preferentially covers water, as both go to low lying areas
      if (pd_(i,0) >= land_cover_.water_transition_depth) {
        res_(i,1) = 1 - res_(i,2);
      } else if (pd_(i,0) <= 0.) {
        res_(i,1) = 0.;
      } else {
        double water_covered = pd_(i,0) / land_cover_.water_transition_depth;
        if (res_(i,2) > water_covered) {
          res_(i,1) = 0;
        } else {
          res_(i,1) = water_covered - res_(i,2);
        }
      }
    }
    res_(i,0) = 1 - res_(i,1) - res_(i,2);

    // if any area is less than eps, give to others
    // if any area fraction is less than eps, give it to the others
    if (res_(i,0) > 0 && res_(i,0) < min_area_) {
      if (res_(i,1) < min_area_) {
        res_(i,2) = 1.;
        res_(i,1) = 0.;
        res_(i,0) = 0.;
      } else {
        res_(i,1) += res_(i,0) * res_(i,1) / (res_(i,1) + res_(i,2));
        res_(i,2) += res_(i,0) * res_(i,2) / (res_(i,1) + res_(i,2));
        res_(i,0) = 0.;
      }
    } else if (res_(i,1) > 0 && res_(i,1) < min_area_) {
      if (res_(i,2) < min_area_) {
        res_(i,0) = 1.;
        res_(i,1) = 0.;
        res_(i,2) = 0.;
      } else {
        res_(i,0) += res_(i,1) * res_(i,0) / (res_(i,0) + res_(i,2));
        res_(i,2) += res_(i,1) * res_(i,2) / (res_(i,0) + res_(i,2));
        res_(i,1) = 0.;
      }
    } else if (res_(i,2) > 0 && res_(i,2) < min_area_) {
      res_(i,0) += res_(i,2) * res_(i,0) / (res_(i,0) + res_(i,1));
      res_(i,1) += res_(i,2) * res_(i,1) / (res_(i,0) + res_(i,1));
      res_(i,2) = 0.;
    }

    assert(fabs(res_(i,0) + res_(i,1) + res_(i,2) - 1.0) < 1.e-10);
    assert(-1.e-10 <= res_(i,0) && res_(i,0) <= 1. + 1.e-10);
    assert(-1.e-10 <= res_(i,1) && res_(i,1) <= 1. + 1.e-10);
    assert(-1.e-10 <= res_(i,2) && res_(i,1) <= 1. + 1.e-10);

    res_(i,0) = fmin(fmax(0., res_(i,0)), 1.);
    res_(i,1) = fmin(fmax(0., res_(i,1)), 1.);
    res_(i,2) = fmin(fmax(0., res_(i,2)), 1.);
  }

  // derivatives not currently provided
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }

 private:
  KeyTag my_key_;
  KeyTag ponded_depth_key_, snow_depth_key_;

  View_type res_;
  cView_type pd_, sd_;

  LandCover land_cover_;
  double min_area_;
};


template <class cView_type, class View_type>
const std::string AreaFractionsThreeComponentModel<cView_type, View_type>::eval_type =
  "area fractions, three components";

template <class cView_type, class View_type>
const std::vector<std::string> AreaFractionsThreeComponentModel<cView_type, View_type>::subfield_names =
  { "ground", "water", "snow" };


using AreaFractionsThreeComponentEvaluator = EvaluatorMultiDOFModelCVByMaterial<AreaFractionsThreeComponentModel>;

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
