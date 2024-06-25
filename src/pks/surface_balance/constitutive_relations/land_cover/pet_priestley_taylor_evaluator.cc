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
#include "seb_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

const std::string PETPriestleyTaylorEvaluator::eval_type =
  "potential evapotranspiration, Priestley-Taylor";


PETPriestleyTaylorEvaluator::PETPriestleyTaylorEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist),
    compatible_(false),
    limiter_(false),
    one_minus_limiter_(false)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  evap_type_ = plist->get<std::string>("evaporation type");
  if (!(evap_type_ == "bare ground" || evap_type_ == "snow" || evap_type_ == "canopy" ||
        evap_type_ == "transpiration")) {
    Errors::Message msg;
    msg << "Priestley-Taylor does not currently support evaporation of type \"" << evap_type_
        << "\".";
    Exceptions::amanzi_throw(msg);
  } else if (evap_type_ == "bare ground") {
    evap_type_ = "ground";
  }

  air_temp_key_ = Keys::readKey(*plist, domain_, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ air_temp_key_, tag });

  surf_temp_key_ = Keys::readKey(*plist, domain_, "surface temperature", "temperature");
  dependencies_.insert(KeyTag{ surf_temp_key_, tag });

  elev_key_ = Keys::readKey(*plist, domain_, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });

  rad_key_ = Keys::readKey(*plist, domain_, "net radiation", "net_radiation");
  dependencies_.insert(KeyTag{ rad_key_, tag });

  limiter_ = plist->get<bool>("include limiter", false);
  if (limiter_) {
    limiter_key_ = Keys::readKey(*plist, domain_, "limiter");
    dependencies_.insert(KeyTag{ limiter_key_, tag });
    limiter_nvecs_ = plist->get<int>("limiter number of dofs", 1);
    limiter_dof_ = plist->get<int>("limiter dof", 0);
  }

  one_minus_limiter_ = plist->get<bool>("include 1 - limiter", false);
  if (one_minus_limiter_) {
    one_minus_limiter_key_ = Keys::readKey(*plist, domain_, "1 - limiter");
    dependencies_.insert(KeyTag{ one_minus_limiter_key_, tag });
    one_minus_limiter_nvecs_ = plist->get<int>("1 - limiter number of dofs", 1);
    one_minus_limiter_dof_ = plist->get<int>("1 - limiter dof", 0);
  }

  land_cover_ = getLandCoverMap(plist->sublist("model parameters"), { "pt_alpha_" + evap_type_ });
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
PETPriestleyTaylorEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  {
    auto air_temp = S.Get<CompositeVector>(air_temp_key_, tag).viewComponent("cell", false);
    auto surf_temp = S.Get<CompositeVector>(surf_temp_key_, tag).viewComponent("cell", false);
    auto elev = S.Get<CompositeVector>(elev_key_, tag).viewComponent("cell", false);
    auto rad = S.Get<CompositeVector>(rad_key_, tag).viewComponent("cell", false);

    auto mesh = result[0]->getMesh();
    auto res = result[0]->viewComponent("cell", false);

    for (const auto& lc : land_cover_) {
      auto lc_ids = mesh->getSetEntities<MemSpace_kind::DEVICE>(
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

      Kokkos::parallel_for(
        "PETPriestleyTaylorEvaluator::Evaluate", lc_ids.extent(0), KOKKOS_LAMBDA(const int i) {
          AmanziMesh::Entity_ID c = lc_ids(i);
          double lh_vap;
          if (is_snow)
            lh_vap = PriestleyTaylor::latentHeatVaporization_snow(air_temp(c, 0));
          else
            lh_vap = PriestleyTaylor::latentHeatVaporization_water(air_temp(c, 0));

          double ps_const = PriestleyTaylor::psychrometricConstant(lh_vap, elev(c, 0));
          double vp_slope = PriestleyTaylor::vaporPressureSlope(air_temp(c, 0));
          double hf_ground = PriestleyTaylor::groundHeatFlux(surf_temp(c, 0), air_temp(c, 0));

          double s1 = vp_slope / (vp_slope + ps_const);
          double s2 = rad(c, 0) - hf_ground; // net radiation balance in W/m^2

          // 1000, density of water kg/m^2/s --> m/s
          res(c, 0) = alpha / lh_vap * s1 * s2 / 1000.;

          // do not allow condensation in P-T
          res(c, 0) = fmax(res(c, 0), 0.0);
        });
    }
  }

  // apply a limiter if requested
  if (limiter_) {
    const auto& limiter = *S.Get<CompositeVector>(limiter_key_, tag).getComponent("cell");
    AMANZI_ASSERT(result.size() == 1);
    AMANZI_ASSERT(result[0]->hasComponent("cell"));
    auto& res = *result[0]->getComponent("cell");
    res.elementWiseMultiply(1, *limiter.getVector(limiter_dof_), res, 0);
  }
  if (one_minus_limiter_) {
    const auto& limiter = *S.Get<CompositeVector>(one_minus_limiter_key_, tag).getComponent("cell");
    AMANZI_ASSERT(result.size() == 1);
    AMANZI_ASSERT(result[0]->hasComponent("cell"));
    auto& res = *result[0]->getComponent("cell");
    res.elementWiseMultiply(-1, *limiter.getVector(one_minus_limiter_dof_), res, 1);
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
    auto limiter = S.Get<CompositeVector>(limiter_key_, tag).viewComponent("cell", false);
    auto evap_val =
      S.Get<CompositeVector>(my_keys_.front().first, tag).viewComponent("cell", false);
    auto res = result[0]->viewComponent("cell", false);
    int limiter_dof(limiter_dof_);

    Kokkos::parallel_for(
      "PETPriestleyTaylorEvaluator::EvaluatePartialDerivative(limiter)",
      res.extent(0),
      KOKKOS_LAMBDA(const int c) {
        double limiter_val = limiter(c, limiter_dof);
        res(c, 0) = limiter_val > 1.e-5 ? evap_val(c, 0) / limiter_val : 0.;
      });

  } else if (one_minus_limiter_ && wrt_key == one_minus_limiter_key_) {
    auto limiter = S.Get<CompositeVector>(one_minus_limiter_key_, tag).viewComponent("cell", false);
    auto evap_val =
      S.Get<CompositeVector>(my_keys_.front().first, tag).viewComponent("cell", false);
    auto res = result[0]->viewComponent("cell", false);
    int one_minus_limiter_dof(one_minus_limiter_dof_);

    Kokkos::parallel_for(
      "PETPriestleyTaylorEvaluator::EvaluatePartialDerivative(limiter)",
      res.extent(0),
      KOKKOS_LAMBDA(const int c) {
        double limiter_val = limiter(c, one_minus_limiter_dof);
        res(c, 0) = limiter_val > 1.e-5 ? -evap_val(c, 0) / (1 - limiter_val) : 0.;
      });
  }
}


void
PETPriestleyTaylorEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!compatible_) {
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
        fac.SetMesh(S.GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
  }
  compatible_ = true;
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
