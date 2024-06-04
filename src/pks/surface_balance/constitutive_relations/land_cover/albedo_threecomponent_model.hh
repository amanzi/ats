/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the two-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A subgrid model for determining albedo for snow and not-snow-covered ground.
/*!

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for three components -- water/ice, land, and
snow.  Note this internally calculates albedo of snow based upon snow density.

Components are: 0 = land, 1 = ice/water, 2 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo-threecomponent-model-spec:
.. admonition:: albedo-threecomponent-model-spec

   * `"albedo ice [-]`" ``[double]`` **0.44**
   * `"albedo water [-]`" ``[double]`` **0.1168**

   * `"emissivity ice [-]`" ``[double]`` **0.98**
   * `"emissivity water [-]`" ``[double]`` **0.995**
   * `"emissivity snow [-]`" ``[double]`` **0.98**

   KEYS:

   - `"subgrid albedos`" **DOMAIN-albedos**
   - `"subgrid emissivities`" **DOMAIN-emissivities**

   DEPENDENCIES:

   - `"snow density`" **SNOW_DOMAIN-density**
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "EvaluatorMultiDOFModelCVByMaterial.hh"
#include "LandCover.hh"
#include "seb_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class AlbedosThreeComponentModel {
 public:
  static const int n_dependencies = 2;
  static const bool provides_derivatives = false;
  static const int n_dofs = 3;
  static const std::vector<std::string> subfield_names;
  static const std::string eval_type;

  explicit AlbedosThreeComponentModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    Key akey = Keys::cleanPListName(*plist);
    auto domain = Keys::getDomain(akey);
    akey = Keys::getVarName(akey);

    Key domain_snow = Keys::readDomainHint(*plist, domain, "surface", "snow");

    Key albedo_key = Keys::in(akey, "albedo") ? akey : "albedos";
    albedo_key_ = Keys::readKeyTag(*plist, domain, "albedos", albedo_key, tag);

    Key emissivity_key = Keys::in(akey, "emissivit") ? akey : "emissivities";
    emissivity_key_ = Keys::readKeyTag(*plist, domain, "emissivities", emissivity_key, tag);

    // dependencies
    unfrozen_fraction_key_ =
      Keys::readKeyTag(*plist, domain, "unfrozen fraction", "unfrozen_fraction", tag);

    // parameters
    is_constant_snow_albedo_ = plist->isParameter("albedo snow [-]");
    if (!is_constant_snow_albedo_) {
      snow_dens_key_ = Keys::readKeyTag(*plist, domain_snow, "snow density", "density", tag);
    }

    a_ice_ = plist->get<double>("albedo ice [-]", 0.44);
    a_water_ = plist->get<double>("albedo water [-]", 0.1168);
    if (is_constant_snow_albedo_) { a_snow_ = plist->get<double>("albedo snow [-]"); }

    e_ice_ = plist->get<double>("emissivity ice [-]", 0.98);
    e_water_ = plist->get<double>("emissivity water [-]", 0.995);
    e_snow_ = plist->get<double>("emissivity snow [-]", 0.98);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    std::string region = model_list.get<std::string>("region");
    land_cover_ = getLandCover(region, model_list, { "albedo_ground", "emissivity_ground" });
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    albedo_ = res[0];
    emissivity_ = res[1];

    uf_ = deps[0];
    if (!is_constant_snow_albedo_) snow_dens_ = deps[1];
  }

  void freeViews()
  {
    albedo_ = View_type();
    emissivity_ = View_type();
    uf_ = cView_type();
    if (!is_constant_snow_albedo_) snow_dens_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { albedo_key_, emissivity_key_ }; }
  KeyTagVector getDependencies() const
  {
    if (is_constant_snow_albedo_) {
      return {
        unfrozen_fraction_key_,
      };
    } else {
      return { unfrozen_fraction_key_, snow_dens_key_ };
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int c) const
  {
    albedo_(c, 0) = land_cover_.albedo_ground;
    albedo_(c, 1) = uf_(c, 0) * a_water_ + (1 - uf_(c, 0)) * a_ice_;
    albedo_(c, 2) = is_constant_snow_albedo_ ? a_snow_ : Functions::albedoSnow(snow_dens_(c, 0));

    emissivity_(c, 0) = land_cover_.emissivity_ground;
    emissivity_(c, 1) = uf_(c, 0) * e_water_ + (1 - uf_(c, 0)) * e_ice_;
    emissivity_(c, 2) = e_snow_;
  }

  // derivatives not currently provided
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }

 private:
  KeyTag albedo_key_, emissivity_key_;
  KeyTag unfrozen_fraction_key_, snow_dens_key_;

  View_type albedo_, emissivity_;
  cView_type uf_, snow_dens_;

  LandCover land_cover_;
  bool is_constant_snow_albedo_;
  double a_ice_, a_water_, a_snow_;
  double e_snow_, e_ice_, e_water_;
};


template <class cView_type, class View_type>
const std::string AlbedosThreeComponentModel<cView_type, View_type>::eval_type =
  "albedos, three components";

template <class cView_type, class View_type>
const std::vector<std::string> AlbedosThreeComponentModel<cView_type, View_type>::subfield_names = {
  "ground",
  "water",
  "snow"
};


using AlbedoThreeComponentEvaluator =
  EvaluatorMultiDOFModelCVByMaterial<AlbedosThreeComponentModel>;

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
