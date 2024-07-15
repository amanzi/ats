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
    Key surf_domain = Keys::readDomainHint(
      *plist, Keys::getDomain(dependencies_.front().first), "domain", "surface");
    surf_krel_key_ = Keys::readKeyTag(
      *plist, surf_domain, "surface relative permeability", "relative_permeability", tag);
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
  // note, we only differentiate the cell quantity here...
  KeyTag wrt{ wrt_key, wrt_tag };
  double rescaling = 1.0 / S.Get<double>("permeability_rescaling", Tags::DEFAULT);

  if (wrt == surf_krel_key_) {
    results[0]->putScalar(0.);
  } else if (wrt == dens_key_) {
    auto res = results[0]->viewComponent("cell", false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent("cell", false);
    auto krel = S.Get<CompositeVector>(krel_key_).viewComponent("cell", false);

    Kokkos::parallel_for(
      "relativepermeabilityevaluator: deriv wrt dens", res.extent(0), KOKKOS_LAMBDA(const int c) {
        res(c, 0) = rescaling * krel(c, 0) / visc(c, 0);
      });

  } else if (wrt == visc_key_) {
    auto res = results[0]->viewComponent("cell", false);
    auto dens = S.Get<CompositeVector>(dens_key_).viewComponent("cell", false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent("cell", false);
    auto krel = S.Get<CompositeVector>(krel_key_).viewComponent("cell", false);

    Kokkos::parallel_for(
      "RelativeHydraulicConductivityEvaluator: deriv wrt visc",
      res.extent(0),
      KOKKOS_LAMBDA(const int c) {
        res(c, 0) = -rescaling * dens(c, 0) * krel(c, 0) / (visc(c, 0) * visc(c, 0));
      });

  } else if (wrt == krel_key_) {
    auto res = results[0]->viewComponent("cell", false);
    auto dens = S.Get<CompositeVector>(dens_key_).viewComponent("cell", false);
    auto visc = S.Get<CompositeVector>(visc_key_).viewComponent("cell", false);

    Kokkos::parallel_for(
      "RelativeHydraulicConductivityEvaluator: deriv wrt krel",
      res.extent(0),
      KOKKOS_LAMBDA(const int c) { res(c, 0) = rescaling * dens(c, 0) / visc(c, 0); });
  }

  if (results[0]->hasComponent("boundary_face"))
    results[0]->getComponent("boundary_face", false)->putScalar(0.);
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
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
