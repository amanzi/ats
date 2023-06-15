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
 Key surface_depth_key_, pressure_head_key_, mask_key_;
 Key sw_domain_name_, pipe_domain_name_;
 
 double manhole_radius_;
 double energ_loss_coeff_; // energy losses coefficient at manhole
 double drain_length_; // drain conduit length

 private:
  static Utils::RegisteredFactory<Evaluator,PipeDrainEvaluator> reg_;

};

} //namespace
} //namespace

#endif
