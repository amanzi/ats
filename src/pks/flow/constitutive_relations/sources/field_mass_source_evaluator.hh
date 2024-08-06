/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le
*/

//! Evaluates tile source mass from tile discharge
/*!

Mass sources into stream/river due to tile discharge

.. _source_mass-tiles-spec:
.. admonition:: source_mass-tiles-spec
*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class FieldMassSourceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit FieldMassSourceEvaluator(Teuchos::ParameterList& plist);
  FieldMassSourceEvaluator(const FieldMassSourceEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new FieldMassSourceEvaluator(*this));
  }

  // virtual void EnsureCompatibility(State& S) override;
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }
  
 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                        const Key& wrt_key,
                                        const Tag& wrt_tag,
                                        const std::vector<CompositeVector*>& result) override{};

 protected:
  Key domain_;
  Key cv_key_;
  Key molar_density_key_;
  Key field_src_key_;
  Teuchos::RCP<Function> QC_curve_;

 private:
  static Utils::RegisteredFactory<Evaluator, FieldMassSourceEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
