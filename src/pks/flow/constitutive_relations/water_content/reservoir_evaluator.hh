/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Incorporates fluxes out of a reservoir/water body.
/*!


.. _field-evaluator-type-reservoir-spec:
.. admonition:: field-evaluator-type-reservoir-spec

   * `"reservoir water body region`" ``[string]`` Name of the surface region
     indicating the reservoir area.

   * `"reservoir outlet region`" ``[string]`` Name of the surface region
     indicating the cell downstream of the reservoir outlet.

   DEPENDENCIES:

   - `"water content`"
   - `"cell volume`"

*/


#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ReservoirModel;

class ReservoirEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  ReservoirEvaluator(Teuchos::ParameterList& plist);
  ReservoirEvaluator(const ReservoirEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ReservoirModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override {
    // pass for now -- not clear if we will support this
  }

  void InitializeFromPlist_();

 protected:
  Key wc_key_;
  Key cv_key_;

  std::string region_waterbody_;
  std::string region_outlet_;

  Teuchos::RCP<ReservoirModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,ReservoirEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

