/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
//! Calculates a face value from cell values through upwinding.

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EvaluatorSecondary.hh"
#include "Factory.hh"

namespace Amanzi {

class EvaluatorUpwindPotentialDifference : public EvaluatorSecondary {
 public:
  explicit EvaluatorUpwindPotentialDifference(Teuchos::ParameterList& plist);

  EvaluatorUpwindPotentialDifference(const EvaluatorUpwindPotentialDifference& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorUpwindPotentialDifference(*this));
  }

  virtual EvaluatorUpwindPotentialDifference& operator=(const Evaluator& other) override;
  EvaluatorUpwindPotentialDifference& operator=(
      const EvaluatorUpwindPotentialDifference& other) = default;


  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key,
                                   const Key& wrt_tag) const override
  {
    // note, this is a bit different than most -- differentiable with respect
    // to the cell quantity only, and if so, we will simply upwind the cell's
    // derivative.
    return S.GetEvaluator(dependencies_[0].first, dependencies_[0].second)
        .IsDifferentiableWRT(S, wrt_key, wrt_tag);
  }

  
  virtual void EnsureCompatibility(State& S) override;

  virtual bool
  UpdateDerivative(State& S, const Key& request, const Key& wrt_key,
                   const Key& wrt_tag) override;

 protected:
  
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

  virtual void Upwind_(const CompositeVector& cells, const CompositeVector& potential,
                       CompositeVector& faces) const;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorUpwindPotentialDifference> reg_;
  
};
  
} // namespace Amanzi
