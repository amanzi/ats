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

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfPumpEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfPumpEvaluator(Teuchos::ParameterList& plist);
  SurfPumpEvaluator(const SurfPumpEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SurfPumpEvaluator(*this));
  }

  // virtual void EnsureCompatibility(State& S) override;
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // calculate of derivatives of this is a tricky thing to do, with
    // non-cell-local terms due to rescaling.  Just turn off derivatives
    // instead.
    return false;
  }

 protected:

  virtual void EnsureCompatibility_Structure_(State& S) override;

  // note, we override Update here because we are working around the fact that
  // this is not really a Monotype evaluator, but also evaluates a flag
  // (pump_on).  Therefore we must override Update to get the non-const flag
  // from State prior to calling Evaluate_.
  virtual void Update_(State& S) override;

  // Two Evaluates are implemented here -- one is the one that accepts the
  // flag, the other is the one required by the EvaluatorSecondaryMonotypeCV
  // API, which is just empty and throws an error.
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result, int& pump_on);

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override {
    AMANZI_ASSERT(false);
  }
  virtual void EvaluatePartialDerivative_(const State& S,
                                        const Key& wrt_key,
                                        const Tag& wrt_tag,
                                        const std::vector<CompositeVector*>& result) override{};

 protected:
  Key cv_key_;
  Key pd_key_;
  Key liq_den_key_;
  Key wc_key_;
  Key pe_key_;

  Key pump_on_key_;

  std::string pump_outlet_region_;
  std::string pump_inlet_region_;
  std::string on_off_region_;
  double max_elev_pumpline_;
  double stage_on_;
  double stage_off_;

  Teuchos::RCP<Function> Q_pump_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfPumpEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
