/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

A Water Retention Model for permafrost requires the computation of saturations
of gas, liquid, and ice.  Multiple models for this are available, based on the
work of `(Painter & Karra 2014) <https://doi.org/10.2136/vzj2013.04.0071>`_

`"evaluator type`" = `"water retention model with ice`"

.. _evaluator-water-retention-model-with-ice-spec:
.. admonition:: evaluator-water-retention-model-with-ice-spec

   ONE OF

   * `"model parameters`" ``[string]`` **"WRM parameters"** Copies this list from
     :ref:`State`'s `"model parameters`" list.

   OR

   * `"model parameters`" ``[wrm-partition-typedinline-spec-list]``

   END

   ONE OF

   * `"permafrost model parameters`" ``[string]`` Copies this list from
     :ref:`State`'s `"model parameters`" list.

   OR

   * `"permafrost model parameters`" ``[permafrost-wrm-partition-typedinline-spec]``

   END


   KEYS:

   - `"liquid saturation`"
   - `"gas saturation`"
   - `"ice saturation`"

   DEPENDENCIES:

   - `"gas-liquid capillary pressure`" **capillary_pressure_gas_liq**
   - `"liquid-ice capillary pressure`" **capillary_pressure_liq_ice**


Like the WRM partition, the permafrost WRM partition provides (region, model)
pairs:

.. _permafrost-wrm-partition-typedinline-spec:
.. admonition:: permafrost-wrm-partition-typedinline-spec

   * `"region`" ``[string]`` region name
   * `"permafrost wrm type`" ``[string]`` name of the model
   * `"_PERMAFROST_WRM_type_ parameters`" ``[_PERMAFROST_WRM_type_-spec]`` See
     below for each type's spec.

*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "wrm_permafrost_model.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMPermafrostEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit WRMPermafrostEvaluator(Teuchos::ParameterList& plist);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMPermafrostModelPartition>& models);
  WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }
  Teuchos::RCP<WRMPermafrostModelPartition> get_WRMPermafrostModels() { return permafrost_models_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    EnsureCompatibility_StructureSame_(S);
  }

  void InitializeFromPlist_();

 protected:
  Key pc_liq_key_;
  Key pc_ice_key_;

  Teuchos::RCP<WRMPermafrostModelPartition> permafrost_models_;
  Teuchos::RCP<WRMPartition> wrms_;

 private:
  static Utils::RegisteredFactory<Evaluator, WRMPermafrostEvaluator> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
