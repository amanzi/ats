/*
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
   - `"porosity`" **DOMAIN_SUB-porosity**
   - `"potential evaporation`" **DOMAIN_SUB-potential_evaporation**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporationDownregulationModel;

class EvaporationDownregulationEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist);
  EvaporationDownregulationEvaluator(const EvaporationDownregulationEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  void InitializeFromPlist_();

  Key sat_gas_key_;
  Key poro_key_;
  Key pot_evap_key_;

  Key domain_surf_;
  Key domain_sub_;

  bool consistent_;

  LandCoverMap land_cover_;
  std::map<std::string, Teuchos::RCP<EvaporationDownregulationModel>> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EvaporationDownregulationEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
