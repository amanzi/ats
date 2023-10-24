/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Downregulates evaporation via vapor diffusion through a dessicated zone.
/*!

Calculates evaporative resistance through a dessicated zone.

Sakagucki and Zeng 2009 equations 9 and 10.

Requires the use of LandCover types, for dessicated zone thickness and Clapp &
Hornberger b.

.. _evaporation-downregulation-evaluator-spec:
.. admonition:: evaporation-downregulation-evaluator-spec

   KEYS:

   - `"saturation gas`" **DOMAIN_SUB-saturation_gas**
   - `"saturation liquid`" **DOMAIN_SUB-saturation_liquid**
   - `"porosity`" **DOMAIN_SUB-porosity**
   - `"potential evaporation`" **DOMAIN_SUB-potential_evaporation**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporationDownregulationModel;

class EvaporationDownregulationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist);
  EvaporationDownregulationEvaluator(const EvaporationDownregulationEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // this will mostly be differentiated with respect to pressure for flow
    // Jacobians, but none of the terms that _really_ depend on p are actually
    // implemented.  That would require differentiating RSoil with respect to
    // s_l, s_g, etc.  But only derivatives wrt potential evaporation are
    // implemented.  That will rarely if ever be p-dependent.  Therefore, this
    // is just turned off to avoid lengthy calculations with 0.
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  void InitializeFromPlist_();

  Key sat_gas_key_;
  Key sat_liq_key_;
  Key poro_key_;
  Key pot_evap_key_;

  Key domain_surf_;
  Key domain_sub_;

  bool consistent_;

  LandCoverMap land_cover_;
  std::map<std::string, Teuchos::RCP<EvaporationDownregulationModel>> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaporationDownregulationEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
