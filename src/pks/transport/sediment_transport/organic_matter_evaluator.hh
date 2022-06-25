/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_ORGANICMATTER_EVALUATOR_
#define AMANZI_ORGANICMATTER_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class OrganicMatterRateEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  OrganicMatterRateEvaluator(Teuchos::ParameterList& plist);

  OrganicMatterRateEvaluator(const OrganicMatterRateEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;
  
  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  //virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

protected:

    // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;


  double Bmax_;
  double Q_db0_;
  Key biomass_key_;

  static Utils::RegisteredFactory<Evaluator,OrganicMatterRateEvaluator> factory_;

};

} //namespace

#endif
