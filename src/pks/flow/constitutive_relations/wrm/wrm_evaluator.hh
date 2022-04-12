/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates saturation through water retention models.
/*!

Water Retention Models (WRMs) determine the saturation as a function of
pressure and the relative permeability as a function of saturation.  Most
commonly used in practice is the van Genuchten model, but others are available.

.. _wrm-evaluator-spec
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
#include "secondary_variables_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  WRMEvaluator(Teuchos::ParameterList& plist);
  WRMEvaluator(Teuchos::ParameterList& plist,
               const Teuchos::RCP<WRMPartition>& wrms);
  WRMEvaluator(const WRMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Teuchos::RCP<WRMPartition> wrms_;
  bool calc_other_sat_;
  Key cap_pres_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> factory_;
  static Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> factory2_;

};

} //namespace
} //namespace

#endif
