/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh.hh"
#include "MeshHelpers_decl.hh"
#include "CompositeVector.hh"
#include "State.hh"

#include "evaluator_upwinding.hh"

namespace Amanzi {
namespace Relations {

EvaluatorUpwinding::EvaluatorUpwinding(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  field_key_ = Keys::readKey(*plist, domain, "field");
  dependencies_.insert(KeyTag{ field_key_, tag });

  strategy_ = plist->get<std::string>("clobber policy", "none");
  if (strategy_ != "none" && strategy_ != "clobber") {
    Errors::Message msg;
    msg << "EvaluatorUpwinding: " << my_keys_.front().first << ": invalid \"clobber policy\" \""
        << strategy_ << "\" -- valid are \"none\" and \"clobber\"";
    Exceptions::amanzi_throw(msg);
  }
}


void
EvaluatorUpwinding::EnsureCompatibility_Structure_(State& S)
{
  for (const auto& keytag : my_keys_) {
    auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(keytag.first, keytag.second);
    my_fac.SetMesh(S.GetMesh(Keys::getDomain(keytag.first)))
      ->SetGhosted()
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  }
}


void
EvaluatorUpwinding::EnsureCompatibility_DerivStructure_(State& S)
{
  for (const auto& keytag : my_keys_) {
    if (S.HasDerivativeSet(keytag.first, keytag.second)) {
      CompositeVectorSpace deriv_fac;
      deriv_fac.SetMesh(S.GetMesh(Keys::getDomain(keytag.first)))
        ->SetGhosted()
        ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 2);

      for (const auto& deriv : S.GetDerivativeSet(keytag.first, keytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first.get());
        S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
          keytag.first, keytag.second, wrt.first, wrt.second, keytag.first)
          .Update(deriv_fac);
      }
    }
  }
}


void
EvaluatorUpwinding::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;

  // the field to upwind lives on cells and (when present) on boundary_face
  // -- boundary faces are then unconditionally seeded from that component
  // before the upwind decision runs (see SeedBoundaryFaces_), since a
  // boundary face has no "other side" cell for the upwind loop to fall
  // back on.
  CompositeVectorSpace field_fac;
  field_fac.SetMesh(S.GetMesh(Keys::getDomain(field_key_)))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(field_key_, tag).Update(field_fac);
}


void
EvaluatorUpwinding::UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag)
{
  std::vector<CompositeVector*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] =
      &S.GetDerivativeW<CompositeVector>(keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    results[j]->putScalarMasterAndGhosted(0.0);
    ++j;
  }

  // The directional dependency (flux_dir/potential) is never differentiated
  // through, so the only dependency whose derivative wrt wrt_key matters
  // here is field_key_.
  Tag tag = my_keys_.front().second;
  if (!S.GetEvaluator(field_key_, tag).IsDifferentiableWRT(S, wrt_key, wrt_tag)) return;

  const auto& dfield_dp = S.GetDerivative<CompositeVector>(field_key_, tag, wrt_key, wrt_tag);

  for (auto* result : results) {
    SeedBoundaryFaces_(S, *result, /*is_derivative=*/true, wrt_key, wrt_tag);
    UpwindDerivative_(S, dfield_dp, *result);
    ApplyClobber_(S, *result, /*is_derivative=*/true, wrt_key, wrt_tag);
  }
}


void
EvaluatorUpwinding::SeedBoundaryFaces_(const State& S,
                                      CompositeVector& result,
                                      bool is_derivative,
                                      const Key& wrt_key,
                                      const Tag& wrt_tag) const
{
  Tag tag = my_keys_.front().second;
  auto mesh = result.getMesh();

  const auto& field = is_derivative ?
    S.GetDerivative<CompositeVector>(field_key_, tag, wrt_key, wrt_tag) :
    S.Get<CompositeVector>(field_key_, tag);
  if (!field.hasComponent("boundary_face")) return;
  const auto& field_bf = *field.getComponent("boundary_face", false);

  if (!is_derivative) {
    const auto& vandelay = mesh->getBoundaryFaceImporter();
    auto& result_f = *result.getComponent("face", false);
    result_f.doExport(field_bf, vandelay, Tpetra::CombineMode::INSERT);
  } else {
    // Writes field_key_'s own "boundary_face" derivative (wrt wrt_key) into
    // slot 0, whatever value State's dependency graph reports for it. Each
    // concrete rule's UpwindDerivative_ reads this value back and places it
    // into the correct final slot (0 for d/d(p_cell), 1 for
    // d/d(p_boundary)) alongside the interior cell's own derivative.
    auto result_f = result.viewComponent("face", false);
    auto field_bf_v = field.viewComponent("boundary_face", false);
    const auto vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);
    const auto face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE, false);
    int nbf = field_bf_v.extent(0);
    for (int bf = 0; bf != nbf; ++bf) {
      int f = face_map->getLocalElement(vandelay_map->getGlobalElement(bf));
      result_f(f, 0) = field_bf_v(bf, 0);
      result_f(f, 1) = 0.0;
    }
  }
}


void
EvaluatorUpwinding::ApplyClobber_(const State& S,
                                  CompositeVector& result,
                                  bool is_derivative,
                                  const Key& wrt_key,
                                  const Tag& wrt_tag) const
{
  if (strategy_ != "clobber") return;
  SeedBoundaryFaces_(S, result, is_derivative, wrt_key, wrt_tag);
}

} // namespace Relations
} // namespace Amanzi
