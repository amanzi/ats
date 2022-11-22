/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation
  2) hydraulic head in the pipe network

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_PIPE_DRAIN_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_PIPE_DRAIN_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PipeDrainEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  explicit
  PipeDrainEvaluator(Teuchos::ParameterList& plist);
  PipeDrainEvaluator(const PipeDrainEvaluator& other) = default;

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
 //What are the names for the flow depth above surface and hydraulic head in pipe network?
 //Key pres_key_, cv_key_;

  double manhole_radius_;
  double energy_losses_coeff_; //at manhole
  double H_max_; //max height of pipe (considering a rectangular cross section)
  double H_; //surface elevation relative to the pipe flow hydraulic head
//we also need some array of cells that says in what cell the manhole is

 private:
  static Utils::RegisteredFactory<Evaluator,PipeDrainEvaluator> reg_;

};

} //namespace
} //namespace

#endif
