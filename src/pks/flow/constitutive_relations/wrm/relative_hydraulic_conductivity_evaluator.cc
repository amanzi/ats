/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "relative_hydraulic_conductivity_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

const std::string RelativeHydraulicConductivityEvaluator::eval_type =
  "relative hydraulic conductivity";

RelativeHydraulicConductivityEvaluator::RelativeHydraulicConductivityEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    use_surface_relperm_(plist->get<bool>("use surface rel perm", false))
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);

  dens_key_ = Keys::readKeyTag(*plist, domain, "density", "molar_density_liquid", tag);
  dependencies_.insert(dens_key_);

  visc_key_ = Keys::readKeyTag(*plist, domain, "viscosity", "viscosity", tag);
  dependencies_.insert(visc_key_);

  krel_key_ =
    Keys::readKeyTag(*plist, domain, "relative permeability", "relative_permeability", tag);
  dependencies_.insert(krel_key_);

  if (use_surface_relperm_) {
    surf_domain_ = Keys::readDomainHint(*plist, domain, "domain", "surface");
    surf_krel_key_ = Keys::readKeyTag(
      *plist, surf_domain_, "surface relative permeability", "relative_permeability", tag);
    dependencies_.insert(surf_krel_key_);
  }
}


Teuchos::RCP<Evaluator>
RelativeHydraulicConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new RelativeHydraulicConductivityEvaluator(*this));
}


void
RelativeHydraulicConductivityEvaluator::Evaluate_(const State& S,
                                                  const std::vector<CompositeVector*>& results)
{
  results[0]->assign(S.Get<CompositeVector>(krel_key_));
  if (use_surface_relperm_) {
    const auto& surf_kr_vec = S.Get<CompositeVector>(surf_krel_key_);
    auto surf_kr = surf_kr_vec.viewComponent("cell", false);
    auto res_bf = results[0]->viewComponent("boundary_face", false);
    const AmanziMesh::MeshCache& mesh = results[0]->getMesh()->getCache();
    const AmanziMesh::MeshCache& surf_mesh = surf_kr_vec.getMesh()->getCache();
    Kokkos::parallel_for(
      "RelativeHydraulicConductivityEvaluator: surf kr to kr",
      surf_kr.extent(0),
      KOKKOS_LAMBDA(const int sc) {
        AmanziMesh::Entity_ID f = surf_mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(mesh, f);
        res_bf(bf, 0) = surf_kr(sc, 0);
      });
  }

  double rescaling = 1.0 / S.Get<double>("permeability_rescaling", Tags::DEFAULT);
  {
    auto dens = S.Get<CompositeVector>(dens_key_).viewComponent("cell", false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent("cell", false);
    auto res = results[0]->viewComponent("cell", false);

    Kokkos::parallel_for(
      "relativepermeabilityevaluator: rho/visc cells", res.extent(0), KOKKOS_LAMBDA(const int c) {
        res(c, 0) = rescaling * dens(c, 0) * res(c, 0) / visc(c, 0);
      });
  }

  if (results[0]->hasComponent("boundary_face")) {
    auto dens = S.Get<CompositeVector>(dens_key_).viewComponent("boundary_face", false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent("boundary_face", false);
    auto res = results[0]->viewComponent("boundary_face", false);

    Kokkos::parallel_for(
      "relativepermeabilityevaluator: rho/visc boundary faces",
      res.extent(0),
      KOKKOS_LAMBDA(const int bf) {
        res(bf, 0) = rescaling * dens(bf, 0) * res(bf, 0) / visc(bf, 0);
      });
  }
}


void
RelativeHydraulicConductivityEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  // coef = rescaling * dens * krel / visc, uniformly on "cell" and
  // "boundary_face" -- which krel value backs "boundary_face" (krel_key_'s
  // own, or surf_krel_key_'s cross-mesh-mapped value) is irrelevant here;
  // it is handled by substituting the dependency derivative that this
  // partial gets multiplied against, in UpdateDerivative_ below.
  KeyTag wrt{ wrt_key, wrt_tag };
  double rescaling = 1.0 / S.Get<double>("permeability_rescaling", Tags::DEFAULT);
  bool have_bf = results[0]->hasComponent("boundary_face");

  for (const auto& comp : std::vector<std::string>{ "cell", "boundary_face" }) {
    if (comp == "boundary_face" && !have_bf) continue;
    auto res = results[0]->viewComponent(comp, false);
    auto dens = S.Get<CompositeVector>(dens_key_).viewComponent(comp, false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent(comp, false);

    if (wrt == dens_key_) {
      auto krel = S.Get<CompositeVector>(krel_key_).viewComponent(comp, false);
      Kokkos::parallel_for(
        "RelativeHydraulicConductivityEvaluator: deriv wrt dens",
        res.extent(0),
        KOKKOS_LAMBDA(const int i) { res(i, 0) = rescaling * krel(i, 0) / visc(i, 0); });

    } else if (wrt == visc_key_) {
      auto krel = S.Get<CompositeVector>(krel_key_).viewComponent(comp, false);
      Kokkos::parallel_for(
        "RelativeHydraulicConductivityEvaluator: deriv wrt visc",
        res.extent(0),
        KOKKOS_LAMBDA(const int i) {
          res(i, 0) = -rescaling * dens(i, 0) * krel(i, 0) / (visc(i, 0) * visc(i, 0));
        });

    } else if (wrt == krel_key_ || wrt == surf_krel_key_) {
      Kokkos::parallel_for(
        "RelativeHydraulicConductivityEvaluator: deriv wrt krel",
        res.extent(0),
        KOKKOS_LAMBDA(const int i) { res(i, 0) = rescaling * dens(i, 0) / visc(i, 0); });
    }
  }
}


void
RelativeHydraulicConductivityEvaluator::UpdateDerivative_(State& S,
                                                          const Key& wrt_key,
                                                          const Tag& wrt_tag)
{
  // Generic chain rule: correct for dens_key_/visc_key_ fully, and for
  // krel_key_'s "cell" contribution. On "boundary_face", it also applies
  // krel_key_'s own derivative (multiplied against the already-correct
  // partial computed above) -- wrong whenever use_surface_relperm_ is true,
  // since Evaluate_ overwrote "boundary_face" with surf_krel_key_'s value
  // instead. Fixed up below.
  EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::UpdateDerivative_(
    S, wrt_key, wrt_tag);

  if (!use_surface_relperm_) return;

  auto& my_deriv =
    S.GetDerivativeW<CompositeVector>(my_keys_.front().first, my_keys_.front().second,
                                      wrt_key, wrt_tag, my_keys_.front().first);
  if (!my_deriv.hasComponent("boundary_face")) return;

  double rescaling = 1.0 / S.Get<double>("permeability_rescaling", Tags::DEFAULT);
  auto res_bf = my_deriv.viewComponent("boundary_face", false);
  auto dens_bf = S.Get<CompositeVector>(dens_key_).viewComponent("boundary_face", false);
  auto visc_bf = S.Get<CompositeVector>(visc_key_).viewComponent("boundary_face", false);
  const AmanziMesh::MeshCache& mesh = my_deriv.getMesh()->getCache();

  Key surf_wrt_key = Keys::getKey(surf_domain_, Keys::getVarName(wrt_key));
  auto& surf_eval = S.GetEvaluator(surf_krel_key_.first, surf_krel_key_.second);

  if (surf_eval.IsDifferentiableWRT(S, surf_wrt_key, wrt_tag)) {
    // substitute surf_krel_key_'s own derivative (cross-mesh-mapped),
    // multiplied against the same ∂coef/∂krel partial already used for
    // the generic "cell" contribution.
    const auto& surf_dkrel =
      S.GetDerivative<CompositeVector>(surf_krel_key_.first, surf_krel_key_.second,
                                       surf_wrt_key, wrt_tag);
    auto surf_dkrel_c = surf_dkrel.viewComponent("cell", false);
    const AmanziMesh::MeshCache& surf_mesh = surf_dkrel.getMesh()->getCache();
    Kokkos::parallel_for(
      "RelativeHydraulicConductivityEvaluator: substitute surf krel deriv (bf)",
      surf_dkrel_c.extent(0),
      KOKKOS_LAMBDA(const int sc) {
        AmanziMesh::Entity_ID f = surf_mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(mesh, f);
        res_bf(bf, 0) = rescaling * dens_bf(bf, 0) * surf_dkrel_c(sc, 0) / visc_bf(bf, 0);
      });
  } else {
    // surf_krel_key_ has zero derivative wrt this wrt_key -- coef_bf has
    // zero dependence on krel at the boundary through this path either.
    Kokkos::parallel_for(
      "RelativeHydraulicConductivityEvaluator: zero krel deriv (bf)",
      res_bf.extent(0),
      KOKKOS_LAMBDA(const int bf) { res_bf(bf, 0) = 0.0; });
  }
}


void
RelativeHydraulicConductivityEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, tag);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.getMesh() != Teuchos::null) {
    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& dep : dependencies_) {
      if (dep != surf_krel_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(my_fac);
      } else {
        CompositeVectorSpace surf_kr_fac;
        surf_kr_fac.SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(surf_kr_fac);
      }
    }
  }

  if (use_surface_relperm_) EnsureCompatibility_SurfDepDerivs_(S);
}


// ---------------------------------------------------------------------------
// surf_krel_key_ depends on a surface-domain key (e.g. surface pressure),
// never the subsurface wrt_key this evaluator's own derivatives are
// requested wrt -- so the generic EnsureCompatibility_DepDerivs_ (which only
// requires a dependency's derivative wrt the SAME wrt_key) never requires
// surf_krel_key_'s derivative. Do so manually here, mapping each wrt_key
// this evaluator has a derivative required for onto its surface-domain
// counterpart.
// ---------------------------------------------------------------------------
void
RelativeHydraulicConductivityEvaluator::EnsureCompatibility_SurfDepDerivs_(State& S)
{
  auto akeytag = my_keys_.front();
  if (!S.HasDerivativeSet(akeytag.first, akeytag.second)) return;

  auto& surf_eval = S.GetEvaluator(surf_krel_key_.first, surf_krel_key_.second);
  for (const auto& deriv : S.GetDerivativeSet(akeytag.first, akeytag.second)) {
    auto wrt = Keys::splitKeyTag(deriv.first.get());
    Key surf_wrt_key = Keys::getKey(surf_domain_, Keys::getVarName(wrt.first));
    if (surf_eval.IsDifferentiableWRT(S, surf_wrt_key, wrt.second)) {
      S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
        surf_krel_key_.first, surf_krel_key_.second, surf_wrt_key, wrt.second);
    }
  }
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
