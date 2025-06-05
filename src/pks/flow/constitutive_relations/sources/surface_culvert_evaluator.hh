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

class SurfCulvertEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfCulvertEvaluator(Teuchos::ParameterList& plist);
  SurfCulvertEvaluator(const SurfCulvertEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SurfCulvertEvaluator(*this));
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
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                        const Key& wrt_key,
                                        const Tag& wrt_tag,
                                        const std::vector<CompositeVector*>& result) override{};

 protected:
  Key domain_;
  Key cv_key_;
  Key pd_key_;
  Key pe_key_;
  Key elev_key_;
  Key liq_den_key_;
  Key wc_key_;

  std::string culvert_inlet_region_;
  std::string culvert_outlet_region_;
  int Nb_;
  double L_;
  double D_;
  double n_;
  double C_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfCulvertEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
