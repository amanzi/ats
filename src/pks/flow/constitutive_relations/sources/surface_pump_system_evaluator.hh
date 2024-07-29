/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! An evaluator for stage-based pump systems 
/*!

This evaluator models stage-based pump station model inside 2D flow area. 
Pump stations can be used to move water between any combination of river reaches, storage areaes or catchment regions. 
Based on pump on-off conditions and pump-operation curve, water is moved from pump-inlet to -outlet region instantly.


.. _surface-pump_system-spec:
.. admonition:: surface-pump_system-spec

   * `"pump inlet region`" ``[str]`` Region of cells from which pump flow is taken out.
   * `"pump outlet region`" ``[str]`` Region of cells from which pump flow is taken out (optional).
   * `"on off reference region`" ``[str]`` Region used to determine when the pump should turn on or off (optional). Defaults to "pump inlet region".
   * `"pump start at stage`" ``[double]`` The water surface elevation that should trigger the pump turning on.
   * `"pump stop at stage`" ``[double]`` The water surface elevation that should shut the pump down.
   * `"maximum pumpline elevation`" ``[double]`` Allows the user to enter the highest elevation in the  pump line. E.g., pumping water over top of a levee.
   * `"pump efficiency  curve`" ``[function]`` This is a function/table of head and flow. Head here is head difference between outlet and inlet (or head by which the water is to be lifted).

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** stage or water surface elevation
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar liquid density`" **DOMAIN-molar_density_liquid** 
   - `"pump on`" **DOMAIN-pump_on_flag** status of pump
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
  Key elev_key_;

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

