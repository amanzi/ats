/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)

*/

//! Simulates water movement through culverts by instant transfer of water from inlet to outlet region.
//! Based on culvert flow equation. Considers, both inlet-controlled and outlet-controlled flow.
/*!
.. _evaluator-surface-culvert-spec:
.. admonition:: evaluator-surface-culvert-spec

   * `"culvert inlet region`" ``[str]`` Region of cells where culvert flow is taken out.
   * `"culvert outlet region`" ``[str]`` Region of cells where culvert flow is introduced.
   * `"number of barrels`" ``[int]`` Number of culvert barrels, default is 1.
   * `"culvert length`" ``[double]`` Length of the culvert in meters, default is 10.
   * `"culvert diameter`" ``[double]`` Diameter of the culvert in meters, default is 1.
   * `"culvert roughness coefficient`" ``[double]`` Manning's roughness coefficient for the culvert, default is 0.013.
   * `"culvert discharge coefficient`" ``[double]`` Discharge coefficient for the culvert, default is 0.6.

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** stage or water surface elevation
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar liquid density`" **DOMAIN-molar_density_liquid** 

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
  double L_feet_;
  double D_feet_;
  double n_;
  double C_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfCulvertEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
