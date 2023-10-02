/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Basic land cover/plant function type
#include "exceptions.hh"
#include "errors.hh"

#include "LandCover.hh"
#include "seb_nan.hh"

namespace Amanzi {
namespace SurfaceBalance {

LandCover::LandCover(Teuchos::ParameterList& plist)
  : rooting_depth_max(plist.get<double>("rooting depth max [m]", NAN)),
    rooting_profile_alpha(plist.get<double>("rooting profile alpha [-]", NAN)),
    rooting_profile_beta(plist.get<double>("rooting profile beta [-]", NAN)),
    stomata_closed_water_potential(
      plist.get<double>("water potential at fully closed stomata [Pa]", NAN)),
    stomata_open_water_potential(
      plist.get<double>("water potential at fully open stomata [Pa]", NAN)),
    leaf_on_doy(plist.get<double>("leaf on time [doy]", NAN)),
    leaf_off_doy(plist.get<double>("leaf off time [doy]", NAN)),
    pt_alpha_snow(plist.get<double>("Priestley-Taylor alpha of snow [-]", NAN)),
    pt_alpha_canopy(plist.get<double>("Priestley-Taylor alpha of canopy [-]", NAN)),
    pt_alpha_ground(plist.get<double>("Priestley-Taylor alpha of bare ground [-]", NAN)),
    pt_alpha_transpiration(plist.get<double>("Priestley-Taylor alpha of transpiration [-]", NAN)),
    albedo_ground(plist.get<double>("albedo of bare ground [-]", NAN)),
    emissivity_ground(plist.get<double>("emissivity of bare ground [-]", NAN)),
    albedo_canopy(plist.get<double>("albedo of canopy [-]", NAN)),
    beers_k_sw(plist.get<double>("Beer's law extinction coefficient, shortwave [-]", NAN)),
    beers_k_lw(plist.get<double>("Beer's law extinction coefficient, longwave [-]", NAN)),
    snow_transition_depth(plist.get<double>("snow transition depth [m]", NAN)),
    water_transition_depth(plist.get<double>("water transition depth [m]", NAN)),
    dessicated_zone_thickness(plist.get<double>("dessicated zone thickness [m]", NAN)),
    clapp_horn_b(plist.get<double>("Clapp and Hornberger b [-]", NAN)),
    rs_method(plist.get<std::string>("soil resistance method", "sakagucki_zeng")),
    roughness_ground(plist.get<double>("roughness length of bare ground [m]", NAN)),
    roughness_snow(plist.get<double>("roughness length of snow [m]", NAN)),
    mannings_n(plist.get<double>("Manning's n [?]", NAN))
{}


LandCoverMap
getLandCover(Teuchos::ParameterList& plist, const std::vector<std::string>& required_pars)
{
  LandCoverMap lcm = Impl::getLandCover(plist);
  for (const auto& lc : lcm) {
    for (const auto& par : required_pars) { Impl::checkValid(lc.first, lc.second, par); }
  }
  return lcm;
}


namespace Impl {

LandCoverMap
getLandCover(Teuchos::ParameterList& plist)
{
  LandCoverMap lc;
  for (auto& item : plist) {
    if (plist.isSublist(item.first)) {
      lc.insert({ item.first, LandCover{ plist.sublist(item.first) } });
    }
  }
  if (lc.size() == 0) {
    Errors::Message message("LandCover is used, but no entries were found in the 'state->initial "
                            "conditions->land cover types' list.");
    Exceptions::amanzi_throw(message);
  }
  return lc;
}

void
throwInvalid(const std::string& region, const std::string& parstr)
{
  Errors::Message msg;
  msg << "LandCover: region \"" << region << "\" missing parameter \"" << parstr << "\"";
  Exceptions::amanzi_throw(msg);
}


void
checkValid(const std::string& region, const LandCover& lc, const std::string& parname)
{
  if (parname == "rooting_depth_max" && std::isnan(lc.rooting_depth_max))
    throwInvalid(region, "rooting depth max [m]");
  if (parname == "rooting_profile_alpha" && std::isnan(lc.rooting_profile_alpha))
    throwInvalid(region, "rooting profile alpha [-]");
  if (parname == "rooting_profile_beta" && std::isnan(lc.rooting_profile_beta))
    throwInvalid(region, "rooting profile beta [-]");

  if (parname == "stomata_closed_water_potential" && std::isnan(lc.stomata_closed_water_potential))
    throwInvalid(region, "water potential at fully closed stomata [Pa]");
  if (parname == "stomata_open_water_potential" && std::isnan(lc.stomata_open_water_potential))
    throwInvalid(region, "water potential at fully open stomata [Pa]");

  if (parname == "pt_alpha_snow" && std::isnan(lc.pt_alpha_snow))
    throwInvalid(region, "Priestley-Taylor alpha of snow [-]");
  if (parname == "pt_alpha_canopy" && std::isnan(lc.pt_alpha_canopy))
    throwInvalid(region, "Priestley-Taylor alpha of canopy [-]");
  if (parname == "pt_alpha_ground" && std::isnan(lc.pt_alpha_ground))
    throwInvalid(region, "Priestley-Taylor alpha of bare ground [-]");
  if (parname == "pt_alpha_transpiration" && std::isnan(lc.pt_alpha_transpiration))
    throwInvalid(region, "Priestley-Taylor alpha of transpiration [-]");

  if (parname == "mannings_n" && std::isnan(lc.mannings_n)) throwInvalid(region, "Manning's n [?]");

  if (parname == "leaf_on_doy" && std::isnan(lc.leaf_on_doy))
    throwInvalid(region, "leaf off time [doy]");
  if (parname == "leaf_off_doy" && std::isnan(lc.leaf_off_doy))
    throwInvalid(region, "leaf off time [doy]");

  if (parname == "emissivity_ground" && std::isnan(lc.emissivity_ground))
    throwInvalid(region, "emissivity of bare ground [-]");
  if (parname == "albedo_ground" && std::isnan(lc.albedo_ground))
    throwInvalid(region, "albedo of bare ground [-]");
  if (parname == "albedo_canopy" && std::isnan(lc.albedo_canopy))
    throwInvalid(region, "albedo of canopy [-]");

  if (parname == "beers_k_lw" && std::isnan(lc.beers_k_lw))
    throwInvalid(region, "Beer's law extinction coefficient, longwave [-]");
  if (parname == "beers_k_sw" && std::isnan(lc.beers_k_sw))
    throwInvalid(region, "Beer's law extinction coefficient, shortwave [-]");

  if (parname == "snow_transition_depth" && std::isnan(lc.snow_transition_depth))
    throwInvalid(region, "snow transition depth [m]");
  if (parname == "water_transition_depth" && std::isnan(lc.water_transition_depth))
    throwInvalid(region, "water transition depth [m]");

  if (parname == "rs_method" && lc.rs_method == std::string("sakagucki_zeng")) {
    if (std::isnan(lc.dessicated_zone_thickness))
      throwInvalid(region, "dessicated zone thickness [m]");
    if (std::isnan(lc.clapp_horn_b)) throwInvalid(region, "Clapp and Hornberger b [-]");
  }

  if (parname == "roughness_ground" && std::isnan(lc.roughness_ground))
    throwInvalid(region, "roughness length of bare ground [m]");
  if (parname == "roughness_snow" && std::isnan(lc.roughness_snow))
    throwInvalid(region, "roughness length of snow [m]");
}

} // namespace Impl
} // namespace SurfaceBalance
} // namespace Amanzi
