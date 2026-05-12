/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Assembles a multicomponent field from per-component sub-fields via column-view aliasing.
/*!

Each component of the assembled field is owned and computed by a dedicated
sub-evaluator whose key is formed by appending ".ComponentName" to the parent
key (e.g. ``tracer2_rate.Tracer2``).  After the parent field's memory is
allocated, the assembler replaces each sub-field's data pointer with a
non-owning column view of the parent.  Writing to a sub-field therefore writes
directly into the correct column of the parent — no copy is needed.

.. _evaluator-component-assembler-spec:
.. admonition:: evaluator-component-assembler-spec

   * `"component names`" ``[Array(string)]`` Names of the sub-components.
     A dependency ``{key}.{name}`` is created automatically for each entry.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ComponentAssemblerEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ComponentAssemblerEvaluator(Teuchos::ParameterList& plist);
  ComponentAssemblerEvaluator(const ComponentAssemblerEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Override Update to wire column views before any sub-evaluator runs.
  bool Update(State& S, const Key& request) override;

 protected:
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
                                  const Key& wrt_key,
                                  const Tag& wrt_tag,
                                  const std::vector<CompositeVector*>& result) override;

  void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  void WireColumnViews_(State& S);

  std::vector<std::string> component_names_;
  bool wired_;

  static Utils::RegisteredFactory<Evaluator, ComponentAssemblerEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi
