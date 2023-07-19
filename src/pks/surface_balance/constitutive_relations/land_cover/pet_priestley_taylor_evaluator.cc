/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

#include "Key.hh"
#include "pet_priestley_taylor_evaluator.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace PriestleyTaylor {

// Convert all units to SI
// input temperature [K], output latent heat [J/kg]
double
latentHeatVaporization_water(double temp_air)
{
  double temp_f = 1.8 * (temp_air - 273.15) + 32; // F
  double lh_vap_calg = 597.3 - (0.5653 * temp_f); // cal/g
  return lh_vap_calg * 4.184 * 1000;              // J/kg
}

double
latentHeatVaporization_snow(double temp_air)
{
  return latentHeatVaporization_water(temp_air); // J/kg
}

// input latent heat [J/kg], elevation [m], output psychrometric const [Pa/C]
double
psychrometricConstant(double lh_vap, double elev)
{
  double elev_ft = elev * 3.281;              // ft
  double lh_vap_calg = lh_vap / 1000 / 4.184; // cal/g
  double psy_kpaC = 1.6286 * (101.3 - (0.003215 * elev_ft)) / lh_vap_calg;
  // Kpa/C
  return psy_kpaC * 1000; // Pa/C
}

// input temperature [K], output slope [Pa/C]
double
vaporPressureSlope(double temp_air)
{
  double tempC = temp_air - 273.15;                            // C
  double vp_sat = Relations::SaturatedVaporPressure(temp_air); // Pa
  return 4098 * vp_sat / std::pow(tempC + 237.3, 2);           // Pa/C
}

// input temperature [K], output heat flux [W/m^2]
double
groundHeatFlux(double temp_ground, double temp_air)
{
  double G = -4.2 * (temp_ground - temp_air); // MJ/m^2/d
  return G * 1e6 / 86400;                     // W/m^2
}

} // namespace PriestleyTaylor


PETPriestleyTaylorEvaluator::PETPriestleyTaylorEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist),
    compatible_(false),
    limiter_(false),
    one_minus_limiter_(false)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  evap_type_ = plist.get<std::string>("evaporation type");
  if (!(evap_type_ == "bare ground" || evap_type_ == "snow" || evap_type_ == "canopy" ||
        evap_type_ == "transpiration")) {
    Errors::Message msg;
    msg << "Priestley-Taylor does not currently support evaporation of type \"" << evap_type_
        << "\".";
    Exceptions::amanzi_throw(msg);
  } else if (evap_type_ == "bare ground") {
    evap_type_ = "ground";
  }

  air_temp_key_ = Keys::readKey(plist, domain_, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ air_temp_key_, tag });

  surf_temp_key_ = Keys::readKey(plist, domain_, "surface temperature", "temperature");
  dependencies_.insert(KeyTag{ surf_temp_key_, tag });

  elev_key_ = Keys::readKey(plist, domain_, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });

  rad_key_ = Keys::readKey(plist, domain_, "net radiation", "net_radiation");
  dependencies_.insert(KeyTag{ rad_key_, tag });

  limiter_ = plist.get<bool>("include limiter", false);
  if (limiter_) {
    limiter_key_ = Keys::readKey(plist, domain_, "limiter");
    dependencies_.insert(KeyTag{ limiter_key_, tag });
    limiter_nvecs_ = plist.get<int>("limiter number of dofs", 1);
    limiter_dof_ = plist.get<int>("limiter dof", 0);
  }

  one_minus_limiter_ = plist.get<bool>("include 1 - limiter", false);
  if (one_minus_limiter_) {
    one_minus_limiter_key_ = Keys::readKey(plist, domain_, "1 - limiter");
    dependencies_.insert(KeyTag{ one_minus_limiter_key_, tag });
    one_minus_limiter_nvecs_ = plist.get<int>("1 - limiter number of dofs", 1);
    one_minus_limiter_dof_ = plist.get<int>("1 - limiter dof", 0);
  }
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
PETPriestleyTaylorEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& air_temp = *S.Get<CompositeVector>(air_temp_key_, tag).ViewComponent("cell", false);
  const auto& surf_temp = *S.Get<CompositeVector>(surf_temp_key_, tag).ViewComponent("cell", false);
  const auto& elev = *S.Get<CompositeVector>(elev_key_, tag).ViewComponent("cell", false);
  const auto& rad = *S.Get<CompositeVector>(rad_key_, tag).ViewComponent("cell", false);

  auto mesh = result[0]->Mesh();
  auto& res = *result[0]->ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    double alpha = 0.;
    bool is_snow = false;
    if (evap_type_ == "snow") {
      alpha = lc.second.pt_alpha_snow;
      is_snow = true;
    } else if (evap_type_ == "ground") {
      alpha = lc.second.pt_alpha_ground;
    } else if (evap_type_ == "canopy") {
      alpha = lc.second.pt_alpha_canopy;
    } else if (evap_type_ == "transpiration") {
      alpha = lc.second.pt_alpha_transpiration;
    }

    for (auto c : lc_ids) {
      double lh_vap;
      if (is_snow)
        lh_vap = PriestleyTaylor::latentHeatVaporization_snow(air_temp[0][c]);
      else
        lh_vap = PriestleyTaylor::latentHeatVaporization_water(air_temp[0][c]);

      double ps_const = PriestleyTaylor::psychrometricConstant(lh_vap, elev[0][c]);
      double vp_slope = PriestleyTaylor::vaporPressureSlope(air_temp[0][c]);
      double hf_ground = PriestleyTaylor::groundHeatFlux(surf_temp[0][c], air_temp[0][c]);

      double s1 = vp_slope / (vp_slope + ps_const);
      double s2 = rad[0][c] - hf_ground; // net radiation balance in W/m^2

      res[0][c] = alpha / lh_vap * s1 * s2 / 1000.; // 1000, density of water
                                                    // kg/m^2/s --> m/s
      // do not allow condensation in P-T
      res[0][c] = std::max(res[0][c], 0.0);
    }
  }

  // apply a limiter if requested
  if (limiter_) {
    const auto& limiter = *S.Get<CompositeVector>(limiter_key_, tag).ViewComponent("cell", false);
#ifdef ENABLE_DBC
    double limiter_max, limiter_min;
    limiter(limiter_dof_)->MaxValue(&limiter_max);
    limiter(limiter_dof_)->MinValue(&limiter_min);
    AMANZI_ASSERT(limiter_max <= 1 + 1e-10);
    AMANZI_ASSERT(limiter_min >= -1e-10);
#endif
    res(0)->Multiply(1, *limiter(limiter_dof_), *res(0), 0);
  }
  if (one_minus_limiter_) {
    const auto& limiter =
      *S.Get<CompositeVector>(one_minus_limiter_key_, tag).ViewComponent("cell", false);
#ifdef ENABLE_DBC
    double limiter_max, limiter_min;
    limiter(one_minus_limiter_dof_)->MaxValue(&limiter_max);
    limiter(one_minus_limiter_dof_)->MinValue(&limiter_min);
    AMANZI_ASSERT(limiter_max <= 1 + 1e-10);
    AMANZI_ASSERT(limiter_min >= -1e-10);
#endif
    res(0)->Multiply(-1, *limiter(one_minus_limiter_dof_), *res(0), 1);
  }
}


void
PETPriestleyTaylorEvaluator::EvaluatePartialDerivative_(const State& S,
                                                        const Key& wrt_key,
                                                        const Tag& wrt_tag,
                                                        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (limiter_ && wrt_key == limiter_key_) {
    const auto& limiter = *S.Get<CompositeVector>(limiter_key_, tag).ViewComponent("cell", false);
    const auto& evap_val =
      *S.Get<CompositeVector>(my_keys_.front().first, tag).ViewComponent("cell", false);
    auto& res_c = *(*result[0]->ViewComponent("cell", false))(0);
    res_c.ReciprocalMultiply(1, *limiter(limiter_dof_), *evap_val(0), 0);
    for (int c = 0; c != res_c.MyLength(); ++c) {
      if (limiter[limiter_dof_][c] < 1.e-5) { res_c[c] = 0.; }
    }
  } else if (one_minus_limiter_ && wrt_key == one_minus_limiter_key_) {
    const auto& limiter =
      *S.Get<CompositeVector>(one_minus_limiter_key_, tag).ViewComponent("cell", false);
    const auto& evap_val =
      *S.Get<CompositeVector>(my_keys_.front().first, tag).ViewComponent("cell", false);
    auto& res_c = *(*result[0]->ViewComponent("cell", false))(0);
    res_c.ReciprocalMultiply(-1, *limiter(one_minus_limiter_dof_), *evap_val(0), 0);
    for (int c = 0; c != res_c.MyLength(); ++c) {
      if (limiter[one_minus_limiter_dof_][c] < 1.e-5) { res_c[c] = 0.; }
    }
  }
}


void
PETPriestleyTaylorEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!compatible_) {
    land_cover_ =
      getLandCover(S.ICList().sublist("land cover types"), { "pt_alpha_" + evap_type_ });

    Tag tag = my_keys_.front().second;
    for (auto& dep : dependencies_) {
      auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, tag);
      if (dep.first == limiter_key_) {
        fac.SetMesh(S.GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, limiter_nvecs_);
      } else if (dep.first == one_minus_limiter_key_) {
        fac.SetMesh(S.GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, one_minus_limiter_nvecs_);
      } else {
        fac.SetMesh(S.GetMesh(domain_))->SetGhosted()->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
  }
  compatible_ = true;
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
