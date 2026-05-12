/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "ComponentAssemblerEvaluator.hh"
#include "errors.hh"
#include "dbc.hh"

namespace Amanzi {
namespace Relations {

ComponentAssemblerEvaluator::ComponentAssemblerEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist),
    wired_(false)
{
  if (!plist_.isParameter("component names")) {
    Errors::Message msg;
    msg << "ComponentAssemblerEvaluator for \"" << my_keys_.front().first
        << "\": missing required parameter \"component names\".";
    Exceptions::amanzi_throw(msg);
  }

  auto names = plist_.get<Teuchos::Array<std::string>>("component names");
  component_names_ = names.toVector();

  if (component_names_.empty()) {
    Errors::Message msg;
    msg << "ComponentAssemblerEvaluator for \"" << my_keys_.front().first
        << "\": \"component names\" must list at least one entry.";
    Exceptions::amanzi_throw(msg);
  }

  // Build the per-component dependency keys: "{my_key}.{component_name}" at my tag.
  // Any dependencies that were already read from the plist (e.g., via "dependencies")
  // are cleared — the assembler's only dependencies are the sub-component fields.
  dependencies_.clear();
  Key my_key = my_keys_.front().first;
  Tag my_tag = my_keys_.front().second;
  for (const auto& comp : component_names_) {
    dependencies_.insert(KeyTag(my_key + "." + comp, my_tag));
  }
}


Teuchos::RCP<Evaluator>
ComponentAssemblerEvaluator::Clone() const
{
  return Teuchos::rcp(new ComponentAssemblerEvaluator(*this));
}


// Wire column views on the first Update call, before sub-evaluators run.
bool
ComponentAssemblerEvaluator::Update(State& S, const Key& request)
{
  if (!wired_) {
    WireColumnViews_(S);
    wired_ = true;
  }
  return EvaluatorSecondaryMonotypeCV::Update(S, request);
}


// Replace each sub-field's data pointer with a non-owning column view of
// the parent field.  After this point, writing to the sub-field writes
// directly into the corresponding column of the parent.
//
// If a sub-evaluator already ran during initialization (before the assembler
// was called), the values it wrote are copied into the parent column before
// the pointer swap, so that data is not lost.
void
ComponentAssemblerEvaluator::WireColumnViews_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag my_tag = my_keys_.front().second;

  auto parent_cv = S.GetPtrW<CompositeVector>(my_key, my_tag, my_key);

  for (int i = 0; i != static_cast<int>(component_names_.size()); ++i) {
    Key subkey = my_key + "." + component_names_[i];

    // GetVector(i) returns a 1-DoF CV whose component data are non-owning
    // views of column i of the parent's EpetraMultiVector.
    auto col_view = parent_cv->GetVector(i);

    // Hold an owning RCP to keep the old allocation alive across SetPtr.
    // Copy whatever the sub-evaluator may have already written into the
    // column view so the parent field is consistent after the pointer swap.
    auto old_sub_cv = S.GetPtr<CompositeVector>(subkey, my_tag);
    col_view->Update(1.0, *old_sub_cv, 0.0);

    S.SetPtr<CompositeVector>(subkey, my_tag, subkey, col_view);
  }
}


// No-op: sub-evaluators write directly through the column views.
void
ComponentAssemblerEvaluator::Evaluate_(const State& S,
                                       const std::vector<CompositeVector*>& result)
{}


void
ComponentAssemblerEvaluator::EvaluatePartialDerivative_(const State& S,
                                                        const Key& wrt_key,
                                                        const Tag& wrt_tag,
                                                        const std::vector<CompositeVector*>& result)
{
  // Derivatives are handled by the sub-evaluators; the assembler has none of its own.
  AMANZI_ASSERT(false);
}


// Tell each sub-field that it must be a 1-DoF field compatible with the
// parent's mesh and component structure.
void
ComponentAssemblerEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag my_tag = my_keys_.front().second;

  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, my_tag);

  if (my_fac.Mesh() == Teuchos::null) {
    // Parent structure not yet known; defer to a later call.
    return;
  }

  for (int i = 0; i != static_cast<int>(component_names_.size()); ++i) {
    Key subkey = my_key + "." + component_names_[i];
    auto& sub_fac = S.Require<CompositeVector, CompositeVectorSpace>(subkey, my_tag);

    // Propagate the parent's mesh and component topology, but with 1 DoF.
    sub_fac.SetMesh(my_fac.Mesh());
    sub_fac.SetGhosted(my_fac.Ghosted());
    for (const auto& comp : my_fac) {
      sub_fac.AddComponent(comp, my_fac.Location(comp), 1);
    }
  }
}

} // namespace Relations
} // namespace Amanzi
