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

double
readPositiveLandCoverParameter(Teuchos::ParameterList& plist, const std::string& name)
{
  double res = plist.get<double>(name, NAN);
  if (res < 0) {
    Errors::Message msg;
    msg << "Invalid land cover parameter \"" << name << "\" in land cover type \"" << plist.name()
        << "\" -- expecting postive value.";
    Exceptions::amanzi_throw(msg);
  }
  return res;
}

double
readNegativeLandCoverParameter(Teuchos::ParameterList& plist, const std::string& name)
{
  double res = plist.get<double>(name, NAN);
  if (res > 0) {
    Errors::Message msg;
    msg << "Invalid land cover parameter \"" << name << "\" in land cover type \"" << plist.name()
        << "\" -- expecting negative value.";
    Exceptions::amanzi_throw(msg);
  }
  return res;
}

double
readZeroOneLandCoverParameter(Teuchos::ParameterList& plist, const std::string& name)
{
  double res = readPositiveLandCoverParameter(plist, name);
  if (res > 1) {
    Errors::Message msg;
    msg << "Invalid land cover parameter \"" << name << "\" in land cover type \"" << plist.name()
        << "\" -- expecting value in range [0,1].";
    Exceptions::amanzi_throw(msg);
  }
  return res;
}


LandCover::LandCover(Teuchos::ParameterList& plist)
  : rooting_depth_max(readPositiveLandCoverParameter(plist, "rooting depth max [m]")),
    rooting_profile_alpha(readPositiveLandCoverParameter(plist, "rooting profile alpha [-]")),
    rooting_profile_beta(plist.get<double>("rooting profile beta [-]", NAN)),
    stomata_closed_capillary_pressure(
      readPositiveLandCoverParameter(plist, "capillary pressure at fully closed stomata [Pa]")),
    stomata_open_capillary_pressure(
      readPositiveLandCoverParameter(plist, "capillary pressure at fully open stomata [Pa]")),
    maximum_xylem_capillary_pressure(
      readPositiveLandCoverParameter(plist, "maximum xylem capillary pressure [Pa]")),
    leaf_on_doy(plist.get<double>("leaf on time [doy]", NAN)),
    leaf_off_doy(plist.get<double>("leaf off time [doy]", NAN)),
    pt_alpha_snow(readPositiveLandCoverParameter(plist, "Priestley-Taylor alpha of snow [-]")),
    pt_alpha_canopy(readPositiveLandCoverParameter(plist, "Priestley-Taylor alpha of canopy [-]")),
    pt_alpha_ground(
      readPositiveLandCoverParameter(plist, "Priestley-Taylor alpha of bare ground [-]")),
    pt_alpha_transpiration(
      readPositiveLandCoverParameter(plist, "Priestley-Taylor alpha of transpiration [-]")),
    albedo_ground(readZeroOneLandCoverParameter(plist, "albedo of bare ground [-]")),
    emissivity_ground(readZeroOneLandCoverParameter(plist, "emissivity of bare ground [-]")),
    albedo_canopy(readZeroOneLandCoverParameter(plist, "albedo of canopy [-]")),
    beers_k_sw(
      readPositiveLandCoverParameter(plist, "Beer's law extinction coefficient, shortwave [-]")),
    beers_k_lw(
      readPositiveLandCoverParameter(plist, "Beer's law extinction coefficient, longwave [-]")),
    snow_transition_depth(readPositiveLandCoverParameter(plist, "snow transition depth [m]")),
    water_transition_depth(readPositiveLandCoverParameter(plist, "water transition depth [m]")),
    roughness_ground(readPositiveLandCoverParameter(plist, "roughness length of bare ground [m]")),
    roughness_snow(readPositiveLandCoverParameter(plist, "roughness length of snow [m]")),
    mannings_n(readPositiveLandCoverParameter(plist, "Manning's n [?]"))
{}


LandCoverMap
getLandCover(Teuchos::ParameterList plist, const std::vector<std::string>& required_pars)
{
  LandCoverMap lcm = Impl::getLandCover(plist);
  for (const auto& lc : lcm) {
    for (const auto& par : required_pars) {
      Impl::checkValid(lc.first, lc.second, par);
    }
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

  if (parname == "stomata_closed_capillary_pressure" &&
      std::isnan(lc.stomata_closed_capillary_pressure))
    throwInvalid(region, "capillary pressure at fully closed stomata [Pa]");
  if (parname == "stomata_open_capillary_pressure" &&
      std::isnan(lc.stomata_open_capillary_pressure))
    throwInvalid(region, "capillary pressure at fully open stomata [Pa]");
  if (parname == "maximum_xylem_capillary_pressure" &&
      std::isnan(lc.maximum_xylem_capillary_pressure))
    throwInvalid(region, "maximum xylem capillary pressure [Pa]");

  if (parname == "pt_alpha_snow" && std::isnan(lc.pt_alpha_snow))
    throwInvalid(region, "Priestley-Taylor alpha of snow [-]");
  if (parname == "pt_alpha_canopy" && std::isnan(lc.pt_alpha_canopy))
    throwInvalid(region, "Priestley-Taylor alpha of canopy [-]");
  if (parname == "pt_alpha_ground" && std::isnan(lc.pt_alpha_ground))
    throwInvalid(region, "Priestley-Taylor alpha of bare ground [-]");
  if (parname == "pt_alpha_transpiration" && std::isnan(lc.pt_alpha_transpiration))
    throwInvalid(region, "Priestley-Taylor alpha of transpiration [-]");

  if (parname == "mannings_n" && std::isnan(lc.mannings_n)
    ) throwInvalid(region, "Manning's n [?]");

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
  if (parname == "roughness_ground" && std::isnan(lc.roughness_ground))
    throwInvalid(region, "roughness length of bare ground [m]");
  if (parname == "roughness_snow" && std::isnan(lc.roughness_snow))
    throwInvalid(region, "roughness length of snow [m]");
}

} // namespace Impl
} // namespace SurfaceBalance
} // namespace Amanzi
