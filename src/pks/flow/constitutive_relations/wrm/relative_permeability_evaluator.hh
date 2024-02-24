/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

This computes the product,

.. math::
   k = \frac{n}{\mu} k_r

which is the scalar portion of the hydraulic conductance (excludes the absolute
permeability).  Relative permeability, k_r, is computed using a water retention
model, which provides generic methods for this quantity based on internal
parameters (e.g. van Genuchten/Mualem, Brooks & Corey, etc.).

While we call this "relative permeability", it is actually the product above,
and might better be called the scalar hydraulic conductance?  The relative
permeability evaluator layers on a few details relative to the relative
permeability model, including the ability to import a surface relative
permeability for incoming infiltration.

type : `"relative permeability`"

.. _relative-permeability-model-spec
.. admonition:: relative-permeability-model-spec

   * `"model parameters`" ``[wrm-spec-list]``

   KEYS
   - `"saturation`" **DOMAIN-saturation_liquid**
   - `"density`" **DOMAIN-molar_density_liquid**
   - `"viscosity`" **DOMAIN-viscosity**

*/

#pragma once

#include "Key.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

template <class RelPermEval_type>
class RelativePermeabilityEvaluator : public RelPermEval_type {
 public:
  RelativePermeabilityEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  RelativePermeabilityEvaluator(const RelativePermeabilityEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key,
          const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  using RelPermEval_type::eval_type;

 protected:
  bool use_surface_relperm_;
  Key surf_krel_key_;

  bool use_dens_visc_factor_;
  Key dens_key_, visc_key_;

  double rescaling_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, RelativePermeabilityEvaluator<RelPermEval_type>> reg_;

};


template <class RelPermEval_type>
RelativePermeabilityEvaluator<RelPermEval_type>::RelativePermeabilityEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist) :
  RelPermEval_type(plist),
  use_surface_relperm_(plist->get<bool>("use surface rel perm", false)),
  use_dens_visc_(plist_.get<bool>("use density on viscosity in rel perm", true)),
  rescaling_(1.0 / plist->get<double>("permeability rescaling", 1.0))
{
  Tag tag = my_keys_.front().second;

  if (use_surface_relperm_) {
    Key surf_domain = Keys::readDomainHint(*plist, Keys::getDomain(dependencies_.front().first), "domain", "surface");
    surf_krel_key_ = Keys::readKey(*plist, surf_domain, "surface relative permeability", "relative_permeability");
    dependencies_.emplace_back(KeyTag(surf_krel_key_, tag));
  }

  if (use_dens_visc_) {
    dens_key_ = Keys::readKey(*plist, domain, "density", "molar_density_liquid");
    dependencies_.emplace_back(KeyTag(dens_key_, tag));

    visc_key_ = Keys::readKey(*plist, domain, "viscosity", "viscosity");
    dependencies_.emplace_back(KeyTag(visc_key_, tag));
  }
}


template <class RelPermEval_type>
Teuchos::RCP<Evaluator>
RelativePermeabilityEvaluator<RelPermEval_type>::Clone() const {
  return Teuchos::rcp(new RelativePermeabilityEvaluator(*this));
}


template <class RelPermEval_type>
void
RelativePermeabilityEvaluator<RelPermEval_type>::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  RelPermEval_type::Evaluate_(S, results);

  Tag tag = my_keys_.front().second;
  if (use_surface_relperm_) {
    const auto& surf_kr_vec = S.Get<CompositeVector>(surf_krel_key_, tag);
    auto surf_kr = surf_kr_vec.viewComponent("cell", false);
    auto res_bf = results[0]->viewComponent("boundary face", false);
    const auto& mesh = *results[0]->getMesh();
    const auto& surf_mesh = *surf_kr_vec.getMesh();
    Kokkos::parallel_for("RelativePermeabilityEvaluator: surf kr to kr",
                         surf_kr.extent_0(),
                         KOKKOS_LAMBDA(const int sc) {
                           AmanziMesh::Entity_ID f = surf_mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
                           AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(mesh, f);
                           res_bf(bf,0) = surf_kr(sc,0);
                         });
  }

  if (use_dens_visc_) {
    auto dens = S.Get<CompositeVector>(dens_key_, tag).viewComponent("cell", false);
    auto visc = S.Get<CompositeVector>(visc_key_, tag).viewComponent("cell", false);
    auto res = results[0]->viewComponent("cell", false);

    Kokkos::parallel_for("relativepermeabilityevaluator: rho/visc cells",
                         res.extent_0(),
                         KOKKOS_LAMBDA(const int c) {
                           res(c,0) = rescaling_ * dens(c,0) * res(c,0) / visc(c,0);
                         });

    if (results[0]->hasComponent("boundary_face")) {
      auto dens = S.Get<CompositeVector>(dens_key_, tag).viewComponent("boundary_face", false);
      auto visc = S.Get<CompositeVector>(visc_key_, tag).viewComponent("boundary_face", false);
      auto res = results[0]->viewComponent("boundary_face", false);

      Kokkos::parallel_for("relativepermeabilityevaluator: rho/visc boundary faces",
                           res.extent_0(),
                           KOKKOS_LAMBDA(const int bf) {
                             res(bf,0) = rescaling_ * dens(bf,0) * res(bf,0) / visc(bf,0);
                           });
    }

  } else if (rescaling_ != 1.0) {
    results[0]->scale(rescaling_);
  }
}


template <class RelPermEval_type>
void
RelativePermeabilityEvaluator<RelPermEval_type>::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key,
        const Tag& wrt_tag,
        const std::vector<CompositeVector*>& results)
{
  if (wrt_key == surf_krel_key_) {
    results[0]->putScalar(0.);
  } else if (wrt_key == dens_key_) {
    RelPermEval_type::Evaluate_(S, results);

    auto visc = S.Get<CompositeVector>(visc_key_, tag).viewComponent("cell", false);
    auto res = results[0]->viewComponent("cell", false);

    Kokkos::parallel_for("RelativePermeabilityEvaluator: deriv d dens",
                         res.extent_0(),
                         KOKKOS_LAMBDA(const int c) {
                           res(c,0) = rescaling_ * res(c,0) / visc(c,0);
                         });
  } else if (wrt_key == visc_key_) {
    RelPermEval_type::Evaluate_(S, results);

    auto dens = S.Get<CompositeVector>(dens_key_, tag).viewComponent("boundary_face", false);
    auto visc = S.Get<CompositeVector>(visc_key_, tag).viewComponent("cell", false);
    auto res = results[0]->viewComponent("cell", false);

    Kokkos::parallel_for("RelativePermeabilityEvaluator: deriv d visc",
                         res.extent_0(),
                         KOKKOS_LAMBDA(const int c) {
                           res(c,0) = -rescaling_ * dens(c,0) * res(c,0) / (visc(c,0) * visc(c,0));
                         });
  } else {
    RelPermEval_type::EvaluatePartialDerivative_(S, wrt_key, wrt_tag, results);

    if (use_dens_visc_) {
      auto dens = S.Get<CompositeVector>(dens_key_, tag).viewComponent("boundary_face", false);
      auto visc = S.Get<CompositeVector>(visc_key_, tag).viewComponent("cell", false);
      auto res = results[0]->viewComponent("cell", false);

      Kokkos::parallel_for("RelativePermeabilityEvaluator: deriv d sat",
                           res.extent_0(),
                           KOKKOS_LAMBDA(const int c) {
                             res(c,0) = rescaling_ * dens(c,0) * res(c,0) / visc(c,0);
                           });
    } else {
      results[0]->scale(rescaling_);
    }
  }

  if (results[0]->hasComponent("boundary_face")) results[0]->getComponent("boundary_face", false)->putScalar(0.);
}


template <class RelPermEval_type>
void
RelativePermeabilityEvaluator<RelPermEval_type>::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, tag);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& dep : dependencies_) {
      if (dep.first != surf_rel_perm_key_) {
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


