/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet @ ornl.gov)
*/

//! SEBThreeComponentEvaluator: evaluates the Surface Energy Balance model on subgrid units.
/*!

Sets up a collection of patches, for portions of the column covered in snow,
ponded water, and vegetated/bare ground.  The surface energy balance on these
area weighted patches are individually calculated then averaged to form the
total quantities.  All down- and up-scaling of relevant quantities are done
through the area weighting, which is calculated by a minimum threshold in snow
and a depression depth/geometry-based approach for water.  All snow is assumed
to first cover water (likely ice), then cover land, as both water and snow
prefer low-lying depressions due to gravity- and wind-driven redistributions,
respectively.


*/

#include "VerboseObject.hh"
#include "seb_threecomponent_evaluator.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

SEBThreeComponentEvaluator::SEBThreeComponentEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), compatible_(false), plist_(plist)
{
  // determine the domain
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);
  Key dtype = Keys::guessDomainType(domain);
  if (dtype == "surface") {
    domain_ = domain;
    domain_ss_ = Keys::readDomainHint(plist_, domain_, "surface", "subsurface");
    domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  } else if (dtype == "domain") {
    domain_ss_ = domain;
    domain_ = Keys::readDomainHint(plist_, domain_ss_, "domain", "surface");
    domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  } else if (dtype == "snow") {
    domain_snow_ = domain;
    domain_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "surface");
    domain_ss_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "subsurface");
  } else {
    domain_snow_ = plist.get<std::string>("snow domain name", domain_snow_);
    domain_ = plist.get<std::string>("surface domain name", domain_);
    domain_ss_ = plist.get<std::string>("subsurface domain name", domain_ss_);
  }
  my_keys_.clear();

  // my keys
  // -- sources
  water_source_key_ = Keys::readKey(plist, domain_, "surface water source", "water_source");
  my_keys_.emplace_back(KeyTag{ water_source_key_, tag });

  energy_source_key_ =
    Keys::readKey(plist, domain_, "surface energy source", "total_energy_source");
  my_keys_.emplace_back(KeyTag{ energy_source_key_, tag });

  ss_water_source_key_ =
    Keys::readKey(plist, domain_ss_, "subsurface water source", "water_source");
  my_keys_.emplace_back(KeyTag{ ss_water_source_key_, tag });

  ss_energy_source_key_ =
    Keys::readKey(plist, domain_ss_, "subsurface energy source", "total_energy_source");
  my_keys_.emplace_back(KeyTag{ ss_energy_source_key_, tag });

  snow_source_key_ = Keys::readKey(plist, domain_snow_, "snow mass source - sink", "source_sink");
  my_keys_.emplace_back(KeyTag{ snow_source_key_, tag });

  new_snow_key_ = Keys::readKey(plist, domain_snow_, "new snow source", "source");
  my_keys_.emplace_back(KeyTag{ new_snow_key_, tag });

  // diagnostics and debugging
  diagnostics_ = plist.get<bool>("save diagnostic data", false);
  if (diagnostics_) {
    // -- diagnostics
    albedo_key_ = Keys::readKey(plist, domain_, "albedo", "albedo");
    my_keys_.emplace_back(KeyTag{ albedo_key_, tag });
    melt_key_ = Keys::readKey(plist, domain_snow_, "snow melt", "melt");
    my_keys_.emplace_back(KeyTag{ melt_key_, tag });
    evap_key_ = Keys::readKey(plist, domain_, "evaporation", "evaporative_flux");
    my_keys_.emplace_back(KeyTag{ evap_key_, tag });
    snow_temp_key_ = Keys::readKey(plist, domain_snow_, "snow temperature", "temperature");
    my_keys_.emplace_back(KeyTag{ snow_temp_key_, tag });
    qE_sh_key_ = Keys::readKey(plist, domain_, "sensible heat flux", "qE_sensible_heat");
    my_keys_.emplace_back(KeyTag{ qE_sh_key_, tag });
    qE_lh_key_ = Keys::readKey(plist, domain_, "latent heat of evaporation", "qE_latent_heat");
    my_keys_.emplace_back(KeyTag{ qE_lh_key_, tag });
    qE_sm_key_ = Keys::readKey(plist, domain_, "latent heat of snowmelt", "qE_snowmelt");
    my_keys_.emplace_back(KeyTag{ qE_sm_key_, tag });
    qE_lw_out_key_ = Keys::readKey(plist, domain_, "outgoing longwave radiation", "qE_lw_out");
    my_keys_.emplace_back(KeyTag{ qE_lw_out_key_, tag });
    qE_cond_key_ = Keys::readKey(plist, domain_, "conducted energy flux", "qE_conducted");
    my_keys_.emplace_back(KeyTag{ qE_cond_key_, tag });
  }

  // dependencies
  // -- met data
  met_sw_key_ =
    Keys::readKey(plist, domain_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(KeyTag{ met_sw_key_, tag });
  met_lw_key_ =
    Keys::readKey(plist, domain_, "incoming longwave radiation", "incoming_longwave_radiation");
  dependencies_.insert(KeyTag{ met_lw_key_, tag });
  met_air_temp_key_ = Keys::readKey(plist, domain_, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ met_air_temp_key_, tag });
  met_vp_air_key_ = Keys::readKey(plist, domain_, "vapor pressure air", "vapor_pressure_air");
  dependencies_.insert(KeyTag{ met_vp_air_key_, tag });
  met_wind_speed_key_ = Keys::readKey(plist, domain_, "wind speed", "wind_speed");
  dependencies_.insert(KeyTag{ met_wind_speed_key_, tag });
  met_prain_key_ = Keys::readKey(plist, domain_, "precipitation rain", "precipitation_rain");
  dependencies_.insert(KeyTag{ met_prain_key_, tag });
  met_psnow_key_ = Keys::readKey(plist, domain_snow_, "precipitation snow", "precipitation");
  dependencies_.insert(KeyTag{ met_psnow_key_, tag });

  // -- snow properties
  snow_depth_key_ = Keys::readKey(plist, domain_snow_, "snow depth", "depth");
  dependencies_.insert(KeyTag{ snow_depth_key_, tag });
  snow_dens_key_ = Keys::readKey(plist, domain_snow_, "snow density", "density");
  dependencies_.insert(KeyTag{ snow_dens_key_, tag });
  snow_death_rate_key_ = Keys::readKey(plist, domain_snow_, "snow death rate", "death_rate");
  dependencies_.insert(KeyTag{ snow_death_rate_key_, tag });

  // -- skin properties
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ mol_dens_key_, tag });
  mass_dens_key_ = Keys::readKey(plist, domain_, "mass density liquid", "mass_density_liquid");
  dependencies_.insert(KeyTag{ mass_dens_key_, tag });
  ponded_depth_key_ = Keys::readKey(plist, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ ponded_depth_key_, tag });
  unfrozen_fraction_key_ = Keys::readKey(plist, domain_, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ unfrozen_fraction_key_, tag });
  sg_albedo_key_ = Keys::readKey(plist, domain_, "albedos", "albedos");
  dependencies_.insert(KeyTag{ sg_albedo_key_, tag });
  sg_emissivity_key_ = Keys::readKey(plist, domain_, "emissivities", "emissivities");
  dependencies_.insert(KeyTag{ sg_emissivity_key_, tag });
  area_frac_key_ = Keys::readKey(plist, domain_, "area fractions", "area_fractions");
  // explicitly excluded to allow snow_death algorithm to work, see #8
  dependencies_.insert(KeyTag{ area_frac_key_, tag });
  surf_temp_key_ = Keys::readKey(plist, domain_, "temperature", "temperature");
  dependencies_.insert(KeyTag{ surf_temp_key_, tag });
  surf_pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(KeyTag{ surf_pres_key_, tag });

  // -- subsurface properties for evaporating bare soil
  sat_gas_key_ = Keys::readKey(plist, domain_ss_, "gas saturation", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });
  sat_liq_key_ = Keys::readKey(plist, domain_ss_, "liquid saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_liq_key_, tag });
  poro_key_ = Keys::readKey(plist, domain_ss_, "porosity", "porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });
  ss_pres_key_ = Keys::readKey(plist, domain_ss_, "subsurface pressure", "pressure");
  dependencies_.insert(KeyTag{ ss_pres_key_, tag });

  // parameters
  min_wind_speed_ = plist.get<double>("minimum wind speed [m s^-1]", 1.0);
  wind_speed_ref_ht_ = plist.get<double>("wind speed reference height [m]", 2.0);
  AMANZI_ASSERT(wind_speed_ref_ht_ > 0.);
}

void
SEBThreeComponentEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  const Relations::ModelParams params(plist_);

  // collect met data
  const auto& qSW_in = *S.Get<CompositeVector>(met_sw_key_, tag).ViewComponent("cell", false);
  const auto& qLW_in = *S.Get<CompositeVector>(met_lw_key_, tag).ViewComponent("cell", false);
  const auto& air_temp =
    *S.Get<CompositeVector>(met_air_temp_key_, tag).ViewComponent("cell", false);
  const auto& vp_air = *S.Get<CompositeVector>(met_vp_air_key_, tag).ViewComponent("cell", false);
  const auto& wind_speed =
    *S.Get<CompositeVector>(met_wind_speed_key_, tag).ViewComponent("cell", false);
  const auto& Prain = *S.Get<CompositeVector>(met_prain_key_, tag).ViewComponent("cell", false);
  const auto& Psnow = *S.Get<CompositeVector>(met_psnow_key_, tag).ViewComponent("cell", false);

  // collect snow properties
  const auto& snow_volumetric_depth =
    *S.Get<CompositeVector>(snow_depth_key_, tag).ViewComponent("cell", false);
  const auto& snow_dens = *S.Get<CompositeVector>(snow_dens_key_, tag).ViewComponent("cell", false);
  const auto& snow_death_rate =
    *S.Get<CompositeVector>(snow_death_rate_key_, tag).ViewComponent("cell", false);

  // collect skin properties
  const auto& mol_dens = *S.Get<CompositeVector>(mol_dens_key_, tag).ViewComponent("cell", false);
  const auto& mass_dens = *S.Get<CompositeVector>(mass_dens_key_, tag).ViewComponent("cell", false);
  const auto& ponded_depth =
    *S.Get<CompositeVector>(ponded_depth_key_, tag).ViewComponent("cell", false);
  const auto& unfrozen_fraction =
    *S.Get<CompositeVector>(unfrozen_fraction_key_, tag).ViewComponent("cell", false);
  const auto& sg_albedo = *S.Get<CompositeVector>(sg_albedo_key_, tag).ViewComponent("cell", false);
  const auto& emissivity =
    *S.Get<CompositeVector>(sg_emissivity_key_, tag).ViewComponent("cell", false);
  const auto& area_fracs =
    *S.Get<CompositeVector>(area_frac_key_, tag).ViewComponent("cell", false);
  const auto& surf_pres = *S.Get<CompositeVector>(surf_pres_key_, tag).ViewComponent("cell", false);
  const auto& surf_temp = *S.Get<CompositeVector>(surf_temp_key_, tag).ViewComponent("cell", false);

  // collect subsurface properties
  const auto& sat_gas = *S.Get<CompositeVector>(sat_gas_key_, tag).ViewComponent("cell", false);
  const auto& sat_liq = *S.Get<CompositeVector>(sat_liq_key_, tag).ViewComponent("cell", false);
  const auto& poro = *S.Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
  const auto& ss_pres = *S.Get<CompositeVector>(ss_pres_key_, tag).ViewComponent("cell", false);

  // collect output vecs
  auto& water_source = *results[0]->ViewComponent("cell", false);
  auto& energy_source = *results[1]->ViewComponent("cell", false);
  auto& ss_water_source = *results[2]->ViewComponent("cell", false);
  auto& ss_energy_source = *results[3]->ViewComponent("cell", false);
  auto& snow_source = *results[4]->ViewComponent("cell", false);
  auto& new_snow = *results[5]->ViewComponent("cell", false);
  water_source.PutScalar(0.);
  energy_source.PutScalar(0.);
  ss_water_source.PutScalar(0.);
  ss_energy_source.PutScalar(0.);
  snow_source.PutScalar(0.);
  new_snow.PutScalar(0.);

  const auto& mesh = *S.GetMesh(domain_);
  const auto& mesh_ss = *S.GetMesh(domain_ss_);

  Epetra_MultiVector *melt_rate(nullptr), *evap_rate(nullptr), *snow_temp(nullptr);
  Epetra_MultiVector *qE_sh(nullptr), *qE_lh(nullptr), *qE_sm(nullptr);
  Epetra_MultiVector *qE_lw_out(nullptr), *qE_cond(nullptr), *albedo(nullptr);
  if (diagnostics_) {
    albedo = results[6]->ViewComponent("cell", false).get();
    albedo->PutScalar(0.);
    melt_rate = results[7]->ViewComponent("cell", false).get();
    melt_rate->PutScalar(0.);
    evap_rate = results[8]->ViewComponent("cell", false).get();
    evap_rate->PutScalar(0.);
    snow_temp = results[9]->ViewComponent("cell", false).get();
    snow_temp->PutScalar(273.15);
    qE_sh = results[10]->ViewComponent("cell", false).get();
    qE_sh->PutScalar(0.);
    qE_lh = results[11]->ViewComponent("cell", false).get();
    qE_lh->PutScalar(0.);
    qE_sm = results[12]->ViewComponent("cell", false).get();
    qE_sm->PutScalar(0.);
    qE_lw_out = results[13]->ViewComponent("cell", false).get();
    qE_lw_out->PutScalar(0.);
    qE_cond = results[14]->ViewComponent("cell", false).get();
    qE_cond->PutScalar(0.);
  }

  unsigned int ncells = water_source.MyLength();
  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh.getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      // get the top cell
      AmanziMesh::Entity_ID subsurf_f = mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
      auto cells = mesh_ss.getFaceCells(subsurf_f, AmanziMesh::Parallel_kind::OWNED);
      AMANZI_ASSERT(cells.size() == 1);

      // met data structure
      Relations::MetData met;
      met.Z_Us = wind_speed_ref_ht_;
      met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      met.QswIn = qSW_in[0][c];
      met.QlwIn = qLW_in[0][c];
      met.air_temp = air_temp[0][c];
      met.vp_air = vp_air[0][c];
      met.Pr = Prain[0][c];

      // bare ground column
      if (area_fracs[0][c] > 0.) {
        Relations::GroundProperties surf;
        surf.temp = surf_temp[0][c];
        surf.pressure = ss_pres[0][cells[0]];
        surf.roughness = lc.second.roughness_ground;
        surf.density_w = mass_dens[0][c];
        surf.dz = lc.second.dessicated_zone_thickness;
        surf.clapp_horn_b = lc.second.clapp_horn_b;
        surf.rs_method = lc.second.rs_method;
        surf.albedo = sg_albedo[0][c];
        surf.emissivity = emissivity[0][c];
        surf.ponded_depth = 0.; // by definition
        surf.porosity = poro[0][cells[0]];
        surf.saturation_gas = sat_gas[0][cells[0]];
        surf.saturation_liq = sat_liq[0][cells[0]];
        surf.unfrozen_fraction = unfrozen_fraction[0][c];
        surf.water_transition_depth = lc.second.water_transition_depth;

        // must ensure that energy is put into melting snow precip, even if it
        // all melts so there is no snow column
        if (area_fracs[2][c] == 0.) {
          met.Ps = Psnow[0][c];
          surf.snow_death_rate = snow_death_rate[0][c]; // m H20 / s
        } else {
          met.Ps = 0.;
          surf.snow_death_rate = 0.;
        }

        // calculate the surface balance
        const Relations::EnergyBalance eb =
          Relations::UpdateEnergyBalanceWithoutSnow(surf, met, params);
        Relations::MassBalance mb = Relations::UpdateMassBalanceWithoutSnow(surf, params, eb);
        Relations::FluxBalance flux = Relations::UpdateFluxesWithoutSnow(surf, met, params, eb, mb);

        // fQe, Me positive is condensation, water flux positive to surface
        water_source[0][c] += area_fracs[0][c] * flux.M_surf;
        energy_source[0][c] += area_fracs[0][c] * flux.E_surf * 1.e-6; // convert to MW/m^2

        double area_to_volume = mesh.getCellVolume(c) / mesh_ss.getCellVolume(cells[0]);
        double ss_water_source_l =
          flux.M_subsurf * area_to_volume * mol_dens[0][c]; // convert from m/m^2/s to mol/m^3/s
        ss_water_source[0][cells[0]] += area_fracs[0][c] * ss_water_source_l;
        double ss_energy_source_l =
          flux.E_subsurf * area_to_volume * 1.e-6; // convert from W/m^2 to MW/m^3
        ss_energy_source[0][cells[0]] += area_fracs[0][c] * ss_energy_source_l;

        snow_source[0][c] += area_fracs[0][c] * flux.M_snow;
        new_snow[0][c] += area_fracs[0][c] * met.Ps;

        if (vo_.os_OK(Teuchos::VERB_EXTREME))
          *vo_.os() << "CELL " << c << " BARE"
                    << ": Ms = " << flux.M_surf << ", Es = " << flux.E_surf * 1.e-6
                    << ", Mss = " << ss_water_source_l << ", Ess = " << ss_energy_source_l
                    << ", Sn = " << flux.M_snow << std::endl;

        // diagnostics
        if (diagnostics_) {
          (*evap_rate)[0][c] -= area_fracs[0][c] * mb.Me;
          (*qE_sh)[0][c] += area_fracs[0][c] * eb.fQh;
          (*qE_lh)[0][c] += area_fracs[0][c] * eb.fQe;
          (*qE_lw_out)[0][c] += area_fracs[0][c] * eb.fQlwOut;
          (*qE_cond)[0][c] += area_fracs[0][c] * eb.fQc;
          (*albedo)[0][c] += area_fracs[0][c] * surf.albedo;

          if (area_fracs[2][c] == 0.) {
            (*qE_sm)[0][c] += area_fracs[0][c] * eb.fQm;
            (*melt_rate)[0][c] += area_fracs[0][c] * mb.Mm;
            (*snow_temp)[0][c] = 273.15;
          }
        }
      }

      // water column
      if (area_fracs[1][c] > 0.) {
        Relations::GroundProperties surf;
        surf.temp = surf_temp[0][c];
        surf.pressure = surf_pres[0][c];
        surf.roughness = lc.second.roughness_ground;
        surf.density_w = mass_dens[0][c];
        surf.dz = lc.second.dessicated_zone_thickness;
        surf.clapp_horn_b = lc.second.clapp_horn_b;
        surf.rs_method = lc.second.rs_method;
        surf.emissivity = emissivity[1][c];
        surf.albedo = sg_albedo[1][c];
        surf.ponded_depth = std::max(lc.second.water_transition_depth, ponded_depth[0][c]);
        surf.porosity = 1.;
        surf.saturation_gas = 0.;
        surf.saturation_liq = sat_liq[0][cells[0]];
        surf.unfrozen_fraction = unfrozen_fraction[0][c];
        surf.water_transition_depth = lc.second.water_transition_depth;

        // must ensure that energy is put into melting snow precip, even if it
        // all melts so there is no snow column
        if (area_fracs[2][c] == 0.) {
          met.Ps = Psnow[0][c];
          surf.snow_death_rate = snow_death_rate[0][c]; // m H20 / s
        } else {
          met.Ps = 0.;
          surf.snow_death_rate = 0.;
        }

        // calculate the surface balance
        const Relations::EnergyBalance eb =
          Relations::UpdateEnergyBalanceWithoutSnow(surf, met, params);
        const Relations::MassBalance mb = Relations::UpdateMassBalanceWithoutSnow(surf, params, eb);
        Relations::FluxBalance flux = Relations::UpdateFluxesWithoutSnow(surf, met, params, eb, mb);

        // fQe, Me positive is condensation, water flux positive to surface
        water_source[0][c] += area_fracs[1][c] * flux.M_surf;
        energy_source[0][c] += area_fracs[1][c] * flux.E_surf * 1.e-6;

        double area_to_volume = mesh.getCellVolume(c) / mesh_ss.getCellVolume(cells[0]);
        double ss_water_source_l =
          flux.M_subsurf * area_to_volume * mol_dens[0][c]; // convert from m/m^2/s to mol/m^3/s
        ss_water_source[0][cells[0]] += area_fracs[1][c] * ss_water_source_l;
        double ss_energy_source_l =
          flux.E_subsurf * area_to_volume * 1.e-6; // convert from W/m^2 to MW/m^3
        ss_energy_source[0][cells[0]] += area_fracs[1][c] * ss_energy_source_l;

        snow_source[0][c] += area_fracs[1][c] * flux.M_snow;
        new_snow[0][c] += area_fracs[1][c] * met.Ps;

        if (vo_.os_OK(Teuchos::VERB_EXTREME))
          *vo_.os() << "CELL " << c << " WATER"
                    << ": Ms = " << flux.M_surf << ", Es = " << flux.E_surf * 1.e-6
                    << ", Mss = " << ss_water_source_l << ", Ess = " << ss_energy_source_l
                    << ", Sn = " << flux.M_snow << std::endl;

        // diagnostics
        if (diagnostics_) {
          (*evap_rate)[0][c] -= area_fracs[1][c] * mb.Me;
          (*qE_sh)[0][c] += area_fracs[1][c] * eb.fQh;
          (*qE_lh)[0][c] += area_fracs[1][c] * eb.fQe;
          (*qE_lw_out)[0][c] += area_fracs[1][c] * eb.fQlwOut;
          (*qE_cond)[0][c] += area_fracs[1][c] * eb.fQc;
          (*albedo)[0][c] += area_fracs[1][c] * surf.albedo;

          if (area_fracs[2][c] == 0.) {
            (*qE_sm)[0][c] += area_fracs[1][c] * eb.fQm;
            (*melt_rate)[0][c] += area_fracs[1][c] * mb.Mm;
            (*snow_temp)[0][c] = 273.15;
          }
        }
      }

      // snow column
      if (area_fracs[2][c] > 0.) {
        Relations::GroundProperties surf;
        surf.temp = surf_temp[0][c];
        surf.pressure = surf_pres[0][c];
        surf.roughness = lc.second.roughness_ground;
        surf.density_w = mass_dens[0][c];
        surf.dz = lc.second.dessicated_zone_thickness;
        surf.clapp_horn_b = lc.second.clapp_horn_b;
        surf.rs_method = lc.second.rs_method; // does not matter
        surf.emissivity = emissivity[2][c];
        surf.albedo = sg_albedo[2][c];
        surf.ponded_depth = 0;    // does not matter
        surf.saturation_gas = 0.; // does not matter
        surf.saturation_liq = sat_liq[0][cells[0]];
        surf.porosity = 1.;                               // does not matter
        surf.unfrozen_fraction = unfrozen_fraction[0][c]; // does not matter
        surf.water_transition_depth = lc.second.water_transition_depth;

        met.Ps = Psnow[0][c] / area_fracs[2][c];

        Relations::SnowProperties snow;
        // take the snow height to be some measure of average thickness -- use
        // volumetric snow depth divided by the area fraction of snow
        snow.height = snow_volumetric_depth[0][c] / area_fracs[2][c];

        // area_fracs may have been set to 1 for snow depth < snow_ground_trans
        // due to min fractional area option in area_fractions evaluator.
        // Decreasing the tol by 1e-6 is about equivalent to a min fractional
        // area of 1e-5 (the default)
        snow.density = snow_dens[0][c];
        snow.albedo = surf.albedo;
        snow.emissivity = surf.emissivity;
        snow.roughness = lc.second.roughness_snow;

        const Relations::EnergyBalance eb =
          Relations::UpdateEnergyBalanceWithSnow(surf, met, params, snow);
        const Relations::MassBalance mb = Relations::UpdateMassBalanceWithSnow(surf, params, eb);
        Relations::FluxBalance flux =
          Relations::UpdateFluxesWithSnow(surf, met, params, snow, eb, mb);

        // fQe, Me positive is condensation, water flux positive to surface.  Subsurf is 0 because of snow
        water_source[0][c] += area_fracs[2][c] * flux.M_surf;
        energy_source[0][c] +=
          area_fracs[2][c] * flux.E_surf * 1.e-6; // convert to MW/m^2 from W/m^2
        snow_source[0][c] += area_fracs[2][c] * flux.M_snow;
        new_snow[0][c] += (met.Ps + std::max(mb.Me, 0.)) * area_fracs[2][c];

        if (vo_.os_OK(Teuchos::VERB_EXTREME))
          *vo_.os() << "CELL " << c << " SNOW"
                    << ": Ms = " << flux.M_surf << ", Es = " << flux.E_surf * 1.e-6
                    << ", Mss = " << 0. << ", Ess = " << 0. << ", Sn = " << flux.M_snow
                    << std::endl;

        // diagnostics
        if (diagnostics_) {
          (*evap_rate)[0][c] -= area_fracs[2][c] * mb.Me;
          (*qE_sh)[0][c] += area_fracs[2][c] * eb.fQh;
          (*qE_lh)[0][c] += area_fracs[2][c] * eb.fQe;
          (*qE_lw_out)[0][c] += area_fracs[2][c] * eb.fQlwOut;
          (*qE_cond)[0][c] += area_fracs[2][c] * eb.fQc;

          (*qE_sm)[0][c] = area_fracs[2][c] * eb.fQm;
          (*melt_rate)[0][c] = area_fracs[2][c] * mb.Mm;
          (*snow_temp)[0][c] = snow.temp;
          (*albedo)[0][c] += area_fracs[2][c] * surf.albedo;
        }
      }
    }
  }

  // debugging
  if (diagnostics_ && vo_.os_OK(Teuchos::VERB_HIGH)) {
    *vo_.os() << "----------------------------------------------------------------" << std::endl
              << "Surface Balance calculation:" << std::endl;
    std::vector<std::string> vnames;
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vnames.push_back("area fractions");
    vecs.push_back(S.GetPtr<CompositeVector>(area_frac_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.clear();
    vecs.clear();
    vnames.push_back("air_temp");
    vecs.push_back(S.GetPtr<CompositeVector>(met_air_temp_key_, tag).ptr());
    vnames.push_back("vp_air");
    vecs.push_back(S.GetPtr<CompositeVector>(met_vp_air_key_, tag).ptr());
    vnames.push_back("precip_rain");
    vecs.push_back(S.GetPtr<CompositeVector>(met_prain_key_, tag).ptr());
    vnames.push_back("precip_snow");
    vecs.push_back(S.GetPtr<CompositeVector>(met_psnow_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("p_ground");
    vecs.push_back(S.GetPtr<CompositeVector>(surf_pres_key_, tag).ptr());
    vnames.push_back("unfrozen_fraction");
    vecs.push_back(S.GetPtr<CompositeVector>(unfrozen_fraction_key_, tag).ptr());
    vnames.push_back("vol_snow_depth");
    vecs.push_back(S.GetPtr<CompositeVector>(snow_depth_key_, tag).ptr());
    vnames.push_back("snow_death");
    vecs.push_back(S.GetPtr<CompositeVector>(snow_death_rate_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("T_ground");
    vecs.push_back(S.GetPtr<CompositeVector>(surf_temp_key_, tag).ptr());
    vnames.push_back("snow_temp");
    vecs.push_back(S.GetPtr<CompositeVector>(snow_temp_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("inc shortwave radiation");
    vecs.push_back(S.GetPtr<CompositeVector>(met_sw_key_, tag).ptr());
    vnames.push_back("inc longwave radiation");
    vecs.push_back(S.GetPtr<CompositeVector>(met_lw_key_, tag).ptr());
    vnames.push_back("inc latent heat");
    vecs.push_back(S.GetPtr<CompositeVector>(qE_lh_key_, tag).ptr());
    vnames.push_back("inc sensible heat");
    vecs.push_back(S.GetPtr<CompositeVector>(qE_sh_key_, tag).ptr());
    vnames.push_back("out longwave radiation");
    vecs.push_back(S.GetPtr<CompositeVector>(qE_lw_out_key_, tag).ptr());
    vnames.push_back("out conducted energy");
    vecs.push_back(S.GetPtr<CompositeVector>(qE_cond_key_, tag).ptr());
    vnames.push_back("out melting energy");
    vecs.push_back(S.GetPtr<CompositeVector>(qE_sm_key_, tag).ptr());

    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("water_src");
    vecs.push_back(S.GetPtr<CompositeVector>(water_source_key_, tag).ptr());
    vnames.push_back("evap flux");
    vecs.push_back(S.GetPtr<CompositeVector>(evap_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    if (vo_.os_OK(Teuchos::VERB_EXTREME))
      *vo_.os() << "CELL 0 Total = "
                << ": Ms = " << water_source[0][0] << ", Es = " << energy_source[0][0]
                << ", Mss = " << ss_water_source[0][99] << ", Ess = " << ss_energy_source[0][99]
                << ", Sn = " << snow_source[0][0] << std::endl;

    vnames.clear();
    vecs.clear();
    vnames.push_back("mass src");
    vnames.push_back("energy src");
    vnames.push_back("snow src");
    vnames.push_back("new snow");
    vecs.push_back(Teuchos::ptr(results[0]));
    vecs.push_back(Teuchos::ptr(results[1]));
    vecs.push_back(Teuchos::ptr(results[4]));
    vecs.push_back(Teuchos::ptr(results[5]));
    db_->WriteVectors(vnames, vecs, true);

    vnames.clear();
    vecs.clear();
    vnames.push_back("sub mass src");
    vnames.push_back("sub energy src");
    vecs.push_back(Teuchos::ptr(results[2]));
    vecs.push_back(Teuchos::ptr(results[3]));
    db_ss_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();
  }
}

void
SEBThreeComponentEvaluator::EvaluatePartialDerivative_(const State& S,
                                                       const Key& wrt_key,
                                                       const Tag& wrt_tag,
                                                       const std::vector<CompositeVector*>& results)
{
  //AMANZI_ASSERT(false);
}

void
SEBThreeComponentEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!compatible_) {
    auto tag = my_keys_.front().second;

    if (db_ == Teuchos::null) {
      db_ = Teuchos::rcp(new Debugger(S.GetMesh(domain_), my_keys_.front().first, plist_));
    }
    if (db_ss_ == Teuchos::null) {
      Teuchos::ParameterList plist(plist_);
      plist.remove("debug cells", false);
      plist.remove("debug faces", false);
      if (plist.isParameter("subsurface debug cells")) {
        plist.set("debug cells", plist.get<Teuchos::Array<int>>("subsurface debug cells"));
      } else {
      }
      if (plist.isParameter("subsurface debug faces"))
        plist.set("debug faces", plist.get<Teuchos::Array<int>>("subsurface debug faces"));
      db_ss_ = Teuchos::rcp(new Debugger(S.GetMesh(domain_ss_), my_keys_.front().first, plist));
    }

    if (land_cover_.size() == 0)
      land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
                                 { "roughness_snow",
                                   "roughness_ground",
                                   "water_transition_depth",
                                   "dessicated_zone_thickness",
                                   "clapp_horn_b",
                                   "rs_method" });

    // use domain name to set the mesh type
    CompositeVectorSpace domain_fac;
    domain_fac.SetMesh(S.GetMesh(domain_))->SetGhosted()->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    CompositeVectorSpace domain_fac_3;
    domain_fac_3.SetMesh(S.GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 3);

    CompositeVectorSpace domain_fac_ss;
    domain_fac_ss.SetMesh(S.GetMesh(domain_ss_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    CompositeVectorSpace domain_fac_snow;
    domain_fac_snow.SetMesh(S.GetMesh(domain_snow_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    for (const auto& dep : dependencies_) {
      auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, tag);
      if (Keys::getDomain(dep.first) == domain_ss_) {
        fac.Update(domain_fac_ss);
      } else if (dep.first == area_frac_key_ || dep.first == sg_albedo_key_ ||
                 dep.first == sg_emissivity_key_) {
        fac.Update(domain_fac_3);
      } else if (Keys::getDomain(dep.first) == domain_snow_) {
        fac.Update(domain_fac_snow);
      } else {
        fac.Update(domain_fac);
      }
    }

    compatible_ = true;
  }
}


void
SEBThreeComponentEvaluator::EnsureCompatibility_Structure_(State& S)
{
  if (!compatible_) {
    // use domain name to set the mesh type
    CompositeVectorSpace domain_fac_owned;
    domain_fac_owned.SetMesh(S.GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    CompositeVectorSpace domain_fac_owned_snow;
    domain_fac_owned_snow.SetMesh(S.GetMesh(domain_snow_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    CompositeVectorSpace domain_fac_owned_ss;
    domain_fac_owned_ss.SetMesh(S.GetMesh(domain_ss_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    for (const auto& key_tag : my_keys_) {
      if (Keys::getDomain(key_tag.first) == domain_) {
        S.Require<CompositeVector, CompositeVectorSpace>(
           key_tag.first, key_tag.second, key_tag.first)
          .Update(domain_fac_owned);
      } else if (Keys::getDomain(key_tag.first) == domain_snow_) {
        S.Require<CompositeVector, CompositeVectorSpace>(
           key_tag.first, key_tag.second, key_tag.first)
          .Update(domain_fac_owned_snow);
      } else if (Keys::getDomain(key_tag.first) == domain_ss_) {
        S.Require<CompositeVector, CompositeVectorSpace>(
           key_tag.first, key_tag.second, key_tag.first)
          .Update(domain_fac_owned_ss);
      } else {
        AMANZI_ASSERT(false);
      }
    }
    // don't flag compatible_ here -- it will be flagged later in EC_ToDeps
  }
}


void
SEBThreeComponentEvaluator::UpdateFieldDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag)
{
  Errors::Message message(
    "SEBTwoComponentEvaluator: cannot differentiate with respect to anything.");
  Exceptions::amanzi_throw(message);
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
