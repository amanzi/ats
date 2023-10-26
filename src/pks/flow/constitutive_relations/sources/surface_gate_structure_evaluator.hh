/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! Evaluates water/solute source which represent effect of distributed subsurface tiles on overland flow
/*!

.. _surface-distributed-tiles-spec:
.. admonition:: surface-distributed-tiles-spec

   * `"number of ditches`" ``[int]`` Number of ditches, corresponding to the number of unique IDs.

   KEYS:
   - `"accumulated source`" **SUBSURFACE_DOMAIN-accumulated_source** Source to the ditch from the tile.
   - `"catchment ID`" **DOMAIN-catchments_id** ID indicating which ditch a given cell drains to.
   - `"catchment fraction`" **DOMAIN-catchments_id** 1/dL, the fraction describing the length scale used in connecting the pipe to the ditch.
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_RELATIONS_SURFDISTTILES_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_SURFDISTTILES_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfGateEvaluator : public EvaluatorSecondary {
 public:
  explicit SurfGateEvaluator(Teuchos::ParameterList& plist);
  SurfGateEvaluator(const SurfGateEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SurfGateEvaluator(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

 protected:
  Key domain_;
  Key cv_key_;
  Key pd_key_;
  Key dependencies_;
  string dp_region_;
  string gate_intake_region_;
  double Q;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfGateEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
#endif
