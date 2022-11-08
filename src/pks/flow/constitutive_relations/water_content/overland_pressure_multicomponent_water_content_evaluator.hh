/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Daniil Svyatsky
*/

#pragma once

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class OverlandPressureMulticomponentWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  explicit
  OverlandPressureMulticomponentWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandPressureMulticomponentWaterContentEvaluator(const OverlandPressureMulticomponentWaterContentEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

 protected:
  Key pres_key_;
  Key cv_key_;
  Key mass_dens_key_, molar_dens_key_;

  bool bar_;  // bar'd variable indicates this is potentially negative for
              // pressures less than atmospheric
  double rollover_;

 private:
  static Utils::RegisteredFactory<Evaluator,OverlandPressureMulticomponentWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
