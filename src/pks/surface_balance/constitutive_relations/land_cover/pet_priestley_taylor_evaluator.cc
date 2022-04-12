/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

#include "Key.hh"
#include "pet_priestley_taylor_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace PriestleyTaylor {

double
latentHeatVaporization_water(double temp_air)
{
  // convert temperature to Fahrenheit
  double temp_f = 1.8 * (temp_air - 273.15) + 32;
  return 597.3 - (0.5653 * temp_f);
}

double
latentHeatVaporization_snow(double temp_air)
{
  return latentHeatVaporization_water(temp_air);
}


double
psychrometricConstant(double lh_vap, double elev)
{
  // convert elevation [m] to elevation [ft]
  double elev_ft = elev * 3.281;
  return 1.6286 * (101.3 - (0.003215 * elev_ft)) / lh_vap;
}

double
vaporPressureSlope(double temp_air)
{
  // temperature conversion from K to C
  double temp_c = (temp_air - 273.15);
  double x = 17.26939*temp_c / (temp_c + 237.3);
  return 4098 * (0.6108 * std::exp(x)) / std::pow(temp_c+237.3,2);
}

double
groundHeatFlux(double temp_ground, double temp_air)
{
  double G = -4.2 * (temp_ground - temp_air);
  return G * 1e6 / 86400; // convert MJ/m^2/d --> W/m^2
}

} // namespace PriestleyTaylor


PETPriestleyTaylorEvaluator::PETPriestleyTaylorEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist),
  compatible_(false),
  limiter_(false),
  one_minus_limiter_(false)
{
  domain_ = Keys::getDomain(my_key_);

  evap_type_ = plist.get<std::string>("evaporation type");
  if (!(evap_type_ == "bare ground" ||
        evap_type_ == "snow" ||
        evap_type_ == "canopy" ||
        evap_type_ == "transpiration")) {
    Errors::Message msg;
    msg << "Priestley-Taylor does not currently support evaporation of type \"" << evap_type_ << "\".";
    Exceptions::amanzi_throw(msg);
  } else if (evap_type_ == "bare ground") {
    evap_type_ = "ground";
  }

  air_temp_key_ = Keys::readKey(plist, domain_, "air temperature", "air_temperature");
  dependencies_.insert(air_temp_key_);

  surf_temp_key_ = Keys::readKey(plist, domain_, "surface temperature", "temperature");
  dependencies_.insert(surf_temp_key_);

  rel_hum_key_ = Keys::readKey(plist, domain_, "relative humidity", "relative_humidity");
  dependencies_.insert(rel_hum_key_);

  elev_key_ = Keys::readKey(plist, domain_, "elevation", "elevation");
  dependencies_.insert(elev_key_);

  rad_key_ = Keys::readKey(plist, domain_, "net radiation", "net_radiation");
  dependencies_.insert(rad_key_);

  limiter_ = plist.get<bool>("include limiter", false);
  if (limiter_) {
    limiter_key_ = Keys::readKey(plist, domain_, "limiter");
    dependencies_.insert(limiter_key_);
    limiter_nvecs_ = plist.get<int>("limiter number of dofs", 1);
    limiter_dof_ = plist.get<int>("limiter dof", 0);
  }

  one_minus_limiter_ = plist.get<bool>("include 1 - limiter", false);
  if (one_minus_limiter_) {
    one_minus_limiter_key_ = Keys::readKey(plist, domain_, "1 - limiter");
    dependencies_.insert(one_minus_limiter_key_);
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void
PETPriestleyTaylorEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& air_temp = *S->GetFieldData(air_temp_key_)->ViewComponent("cell", false);
  const auto& surf_temp = *S->GetFieldData(surf_temp_key_)->ViewComponent("cell", false);
  const auto& rel_hum = *S->GetFieldData(rel_hum_key_)->ViewComponent("cell", false);
  const auto& elev = *S->GetFieldData(elev_key_)->ViewComponent("cell", false);
  const auto& rad = *S->GetFieldData(rad_key_)->ViewComponent("cell",false);

  auto mesh = result->Mesh();
  auto& res = *result->ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

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
      double lh_vap_si = lh_vap * 4.184 * 1000.; // converts cal/gm to J/kg

      double vp_slope = PriestleyTaylor::vaporPressureSlope(air_temp[0][c]);
      double hf_ground = PriestleyTaylor::groundHeatFlux(surf_temp[0][c], air_temp[0][c]);

      double s1 = vp_slope / (vp_slope + ps_const);
      double s2 = rad[0][c] - hf_ground; // net radiation balance in W/m^2

      res[0][c] = alpha / lh_vap_si * s1 * s2 / 1000.;  // 1000, density of
                                                       // water converts from
                                                       // kg/m^2/s --> m/s
      // do not allow condensation in P-T
      res[0][c] = std::max(res[0][c],0.0);
    }
  }

  // apply a limiter if requested
  if (limiter_) {
    const auto& limiter = *S->GetFieldData(limiter_key_)->ViewComponent("cell", false);
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
    const auto& limiter = *S->GetFieldData(one_minus_limiter_key_)->ViewComponent("cell", false);
#ifdef ENABLE_DBC
    double limiter_max, limiter_min;
    limiter.MaxValue(&limiter_max);
    limiter.MinValue(&limiter_min);
    AMANZI_ASSERT(limiter_max <= 1 + 1e-10);
    AMANZI_ASSERT(limiter_min >= -1e-10);
#endif
    res.Multiply(-1, limiter, res, 1);
  }
}


void
PETPriestleyTaylorEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  if (limiter_ && wrt_key == limiter_key_) {
    const auto& limiter = *S->GetFieldData(limiter_key_)->ViewComponent("cell", false);
    const auto& evap_val = *S->GetFieldData(my_key_)->ViewComponent("cell", false);
    auto& res_c = *(*result->ViewComponent("cell", false))(0);
    res_c.ReciprocalMultiply(1, *limiter(limiter_dof_), *evap_val(0), 0);
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (limiter[limiter_dof_][c] < 1.e-5) {
        res_c[c] = 0.;
      }
    }
  } else if (one_minus_limiter_ && wrt_key == one_minus_limiter_key_) {
    const auto& limiter = *S->GetFieldData(one_minus_limiter_key_)->ViewComponent("cell", false);
    const auto& evap_val = *S->GetFieldData(my_key_)->ViewComponent("cell", false);
    auto& res_c = *result->ViewComponent("cell", false);
    res_c.ReciprocalMultiply(-1, limiter, evap_val, 0);
    for (int c=0; c!=res_c.MyLength(); ++c) {
      if (limiter[0][c] < 1.e-5) {
        res_c[0][c] = 0.;
      }
    }
  }
}


void
PETPriestleyTaylorEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (!compatible_) {
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"pt_alpha_"+evap_type_});

    // see if we can find a master fac
    auto my_fac = S->RequireField(my_key_, my_key_);
    my_fac->SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>("visualize", true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    for (auto dep_key : dependencies_) {
      auto fac = S->RequireField(dep_key);
      if (dep_key == limiter_key_) {
        fac->SetMesh(S->GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, limiter_nvecs_);
      } else {
        fac->SetMesh(S->GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
      }

      // Recurse into the tree to propagate info to leaves.
      S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
    }
  }
  compatible_ = true;
}



} //namespace
} //namespace
} //namespace

