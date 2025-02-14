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

#include "EvaluatorSecondary.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PipeDrainEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit PipeDrainEvaluator(Teuchos::ParameterList& plist);
  PipeDrainEvaluator(const PipeDrainEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;
  virtual void EnsureCompatibility_ToDeps_(State& S, const CompositeVectorSpace& fac) override;
  void CreateCellMap(const State& S);

 protected:
  Key surface_depth_key_, pressure_head_key_, mask_key_;
  Key sw_domain_name_, pipe_domain_name_;
  Key manhole_map_key_;
  Key surface_bathymetry_key_, pipe_bathymetry_key_;

  double manhole_radius_;
  // energy losses coefficients at manhole
  // see Table 3 in "Experimental calibration and validation of sewer/surface
  // flow exchange equations in steady and unsteady flow conditions" by Rubinato et al.
  double energ_loss_coeff_weir_;    // weir
  double energ_loss_coeff_subweir_; // submerged weir
  double energ_loss_coeff_orifice_; // orifice
  double
    sink_source_coeff_; // coefficient that determines sink or source (when using same evaluator file for pipe or surface flow)
  bool pipe_flag_, sw_flag_, pipe_map_created_ = false, sw_map_created_ = false;

  std::vector<int> pipe_map_, sw_map_;

 private:
  static Utils::RegisteredFactory<Evaluator, PipeDrainEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
