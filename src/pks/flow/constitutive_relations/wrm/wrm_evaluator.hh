/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates saturation through water retention models.
/*!

Water Retention Models (WRMs) determine the saturation as a function of
pressure and the relative permeability as a function of saturation.  Most
commonly used in practice is the van Genuchten model, but others are available default default;

`"evaluator type`" = `"wrm`"

.. _wrm-evaluator-spec:
.. admonition:: wrm-evaluator-spec

   * `"WRM parameters`" ``[WRM-typedinline-spec-list]``

   KEYS:

   - `"saturation`" **determined from evaluator name** The name
     of the liquid saturation -- typically this is determined from
     the evaluator name and need not be set.
   - `"other saturation`"  **determined from evaluator name**
     The name of the other saturation, usually gas -- typically this is determined
     from the evaluator name and need not be set.
   - `"capillary pressure`"` **DOMAIN-capillary_pressure_gas_liq**
     The name of the capillary pressure.

*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_

#include "wrm_partition.hh"
#include "wrm.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit WRMEvaluator(Teuchos::ParameterList& plist);
  WRMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  WRMEvaluator(const WRMEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  void InitializeFromPlist_();

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

 protected:
  Teuchos::RCP<WRMPartition> wrms_;
  bool calc_other_sat_;
  Key cap_pres_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WRMEvaluator> factory_;
  static Utils::RegisteredFactory<Evaluator, WRMEvaluator> factory2_;
};

} // namespace Flow
} // namespace Amanzi

#endif
