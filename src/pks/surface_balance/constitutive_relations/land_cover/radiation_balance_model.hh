/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

// See radiation_balance_evaluator.hh for documentation.

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "LandCover.hh"
#include "seb_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class RadiationBalanceModel {
 public:
  static const int n_dependencies = 9;
  static const bool provides_derivatives = false;
  static const std::string eval_type;

  explicit RadiationBalanceModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Key akey = Keys::cleanPListName(*plist);
    Key domain = Keys::getDomain(akey);
    Tag tag(plist->get<std::string>("tag"));

    Key dtype = Keys::guessDomainType(domain);
    Key domain_surf, domain_snow, domain_canopy;
    if (dtype == "surface") {
      domain_surf = domain;
      domain_snow = Keys::readDomainHint(*plist, domain_surf, "surface", "snow");
      domain_canopy = Keys::readDomainHint(*plist, domain_surf, "surface", "canopy");
    } else if (dtype == "canopy") {
      domain_canopy = domain;
      domain_snow = Keys::readDomainHint(*plist, domain_canopy, "canopy", "snow");
      domain_surf = Keys::readDomainHint(*plist, domain_canopy, "canopy", "surface");
    } else if (dtype == "snow") {
      domain_snow = domain;
      domain_canopy = Keys::readDomainHint(*plist, domain_snow, "snow", "canopy");
      domain_surf = Keys::readDomainHint(*plist, domain_snow, "snow", "surface");
    } else {
      domain_surf = plist->get<std::string>("surface domain name");
      domain_snow = plist->get<std::string>("snow domain name");
      domain_canopy = plist->get<std::string>("canopy domain name");
    }
    akey = Keys::getVarName(akey);

    // my keys
    rad_bal_surf_key_ = Keys::readKeyTag(*plist,
                                         domain_surf,
                                         "surface radiation balance",
                                         dtype == "surface" ? akey : "radiation_balance");
    rad_bal_snow_key_ = Keys::readKeyTag(
      *plist, domain_snow, "snow radiation balance", dtype == "snow" ? akey : "radiation_balance");
    rad_bal_canopy_key_ = Keys::readKeyTag(*plist,
                                           domain_canopy,
                                           "canopy radiation balance",
                                           dtype == "canopy" ? akey : "radiation_balance");

    // dependencies
    albedo_surf_key_ = Keys::readKeyTag(*plist, domain_surf, "surface albedos", "albedos", tag);
    emissivity_surf_key_ =
      Keys::readKeyTag(*plist, domain_surf, "surface emissivities", "emissivities", tag);
    sw_in_key_ = Keys::readKeyTag(
      *plist, domain_surf, "incoming shortwave radiation", "incoming_shortwave_radiation", tag);
    lw_in_key_ = Keys::readKeyTag(
      *plist, domain_surf, "incoming longwave radiation", "incoming_longwave_radiation", tag);
    temp_surf_key_ =
      Keys::readKeyTag(*plist, domain_surf, "surface temperature", "temperature", tag);
    temp_snow_key_ = Keys::readKeyTag(*plist, domain_snow, "snow temperature", "temperature", tag);
    temp_canopy_key_ =
      Keys::readKeyTag(*plist, domain_canopy, "canopy temperature", "temperature", tag);
    area_frac_key_ = Keys::readKeyTag(*plist, domain_surf, "area fractions", "area_fractions", tag);
    lai_key_ = Keys::readKeyTag(*plist, domain_canopy, "leaf area index", "leaf_area_index", tag);

    // land cover struct
    std::string region_name = plist->sublist("model parameters").get<std::string>("region");
    lc_ = getLandCover(region_name,
                       plist->sublist("model parameters"),
                       { "albedo_canopy", "beers_k_lw", "beers_k_sw" });
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(res.size() == 3);
    rad_bal_surf = res[0];
    rad_bal_snow = res[1];
    rad_bal_canopy = res[2];

    AMANZI_ASSERT(deps.size() == 9);
    albedo_surf = deps[0];
    emissivity_surf = deps[1];
    sw_in = deps[2];
    lw_in = deps[3];
    temp_surf = deps[4];
    temp_snow = deps[5];
    temp_canopy = deps[6];
    area_frac = deps[7];
    lai = deps[8];
  }

  void freeViews()
  {
    rad_bal_surf = View_type();
    rad_bal_snow = View_type();
    rad_bal_canopy = View_type();
    albedo_surf = cView_type();
    emissivity_surf = cView_type();
    sw_in = cView_type();
    lw_in = cView_type();
    temp_surf = cView_type();
    temp_snow = cView_type();
    temp_canopy = cView_type();
    area_frac = cView_type();
    lai = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return { rad_bal_surf_key_, rad_bal_snow_key_, rad_bal_canopy_key_ };
  }
  KeyTagVector getDependencies() const
  {
    return { albedo_surf_key_, emissivity_surf_key_, sw_in_key_,     lw_in_key_, temp_surf_key_,
             temp_snow_key_,   temp_canopy_key_,     area_frac_key_, lai_key_ };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    // NOTE: emissivity = absorptivity, we use e to notate both
    // Beer's law to find absorptivity of canopy
    double e_can_sw = Functions::beersLawAbsorptivity(lc_.beers_k_sw, lai(i, 0));
    double e_can_lw = Functions::beersLawAbsorptivity(lc_.beers_k_lw, lai(i, 0));

    // sw atm to canopy and surface
    double sw_atm_can = e_can_sw * sw_in(i, 0);
    double sw_atm_surf = sw_in(i, 0) - sw_atm_can;

    // reflected sw (can be important off of snow)
    double sw_grnd_can =
      e_can_sw * sw_atm_surf *
      (albedo_surf(i, 0) * area_frac(i, 0) + albedo_surf(i, 1) * area_frac(i, 1));

    // lw atm to canopy and surface
    double lw_atm_can = e_can_lw * lw_in(i, 0);
    double lw_atm_surf = lw_in(i, 0) - lw_atm_can;

    // lw out of each layer
    double lw_surf = Functions::outgoingLongwaveRadiation(temp_surf(i, 0), emissivity_surf(i, 0));
    double lw_snow = Functions::outgoingLongwaveRadiation(temp_snow(i, 0), emissivity_surf(i, 1));
    double lw_can = Functions::outgoingLongwaveRadiation(temp_canopy(i, 0), e_can_lw);

    // surface connections
    double lw_down = lw_atm_surf + lw_can;
    double lw_up_surf = (1 - emissivity_surf(i, 0)) * lw_down + lw_surf;
    double lw_up_snow = (1 - emissivity_surf(i, 1)) * lw_down + lw_snow;

    // radiation balances -- see Figure 4.1 in CLM Tech Note
    rad_bal_surf(i, 0) = (1 - albedo_surf(i, 0)) * sw_atm_surf + lw_down - lw_up_surf;
    rad_bal_snow(i, 0) = (1 - albedo_surf(i, 1)) * sw_atm_surf + lw_down - lw_up_snow;

    rad_bal_canopy(i, 0) = (1 - lc_.albedo_canopy) * (sw_atm_can + sw_grnd_can) + lw_atm_can -
                           2 * lw_can + area_frac(i, 0) * e_can_lw * lw_up_surf +
                           area_frac(i, 1) * e_can_lw * lw_up_snow;
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<4>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<5>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<6>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<7>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<8>, const int i) const { assert(false); }

 private:
  View_type rad_bal_surf, rad_bal_snow, rad_bal_canopy;
  cView_type albedo_surf, emissivity_surf, sw_in, lw_in;
  cView_type temp_surf, temp_snow, temp_canopy, area_frac, lai;

  KeyTag rad_bal_surf_key_, rad_bal_snow_key_, rad_bal_canopy_key_;
  KeyTag albedo_surf_key_, emissivity_surf_key_, sw_in_key_, lw_in_key_;
  KeyTag temp_surf_key_, temp_snow_key_, temp_canopy_key_, area_frac_key_, lai_key_;

  LandCover lc_;
};


template <class cView_type, class View_type>
const std::string RadiationBalanceModel<cView_type, View_type>::eval_type = "radiation balance";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
