/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Three phase water content: vapor, liquid, and ice.

/*!

.. math::
  Theta = (n_l * s_l + n_i * s_i + n_g * s_g * \omega_g ) * \phi * |E|

* `"porosity key`" ``[string]`` **DOMAIN-porosity** 

* `"molar density liquid key`" ``[string]`` **DOMAIN-molar_density_liquid** 
* `"saturation liquid key`" ``[string]`` **DOMAIN-saturation_liquid** 

* `"molar density ice key`" ``[string]`` **DOMAIN-molar_density_ice** 
* `"saturation ice key`" ``[string]`` **DOMAIN-saturation_ice** 

* `"molar density gas key`" ``[string]`` **DOMAIN-molar_density_gas** 
* `"saturation gas key`" ``[string]`` **DOMAIN-saturation_gas** 
* `"mol frac gas key`" ``[string]`` **DOMAIN-mol_frac_gas** The molar fraction of water vapor in the gaseous phase.

* `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**

EVALUATORS:
- `"porosity`"
- `"molar density liquid`"
- `"saturation liquid`"
- `"molar density ice`"
- `"saturation ice`"
- `"molar density gas`"
- `"saturation gas`"
- `"molar fraction gas`"
- `"cell volume`"

*/
/*


  Generated via evaluator_generator with:
Water content for a three-phase, gas+liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ThreePhaseWaterContentModel;

class ThreePhaseWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  ThreePhaseWaterContentEvaluator(Teuchos::ParameterList& plist);
  ThreePhaseWaterContentEvaluator(const ThreePhaseWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<ThreePhaseWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<ThreePhaseWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ThreePhaseWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
