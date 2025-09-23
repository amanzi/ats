/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Water Retention Models (WRMs) determine the saturation as a function of
pressure and the relative permeability as a function of saturation.  Most
commonly used in practice is the van Genuchten model, but others are available default default;

`"evaluator type`" = `"wrm`"

.. _evaluator-wrm-spec:
.. admonition:: evaluator-wrm-spec

   * `"model parameters`" ``[string]`` **"WRM parameters"** This will
     copy `"WRM parameters`" given in `"model parameters`" under state here to
     evaluate WRM.

   KEYS:

   - `"saturation`" **determined from evaluator name** The name
     of the liquid saturation -- typically this is determined from
     the evaluator name and need not be set.
   - `"other saturation`"  **determined from evaluator name**
     The name of the other saturation, usually gas -- typically this is determined
     from the evaluator name and need not be set.

   DEPENDENCIES:

   - `"capillary pressure`"` **DOMAIN-capillary_pressure_gas_liq**
     The name of the capillary pressure.

The resulting `"model parameters`" are expected to be a WRM partition.  A WRM
partition is a list of (region, WRM) pairs, where the regions partition the
mesh.

.. _wrm-partition-typedinline-spec:
.. admonition:: wrm-partition-typedinline-spec

   * `"region`" ``[string]`` Region on which the WRM is valid.
   * `"WRM type`" ``[string]`` Name of the WRM type.
   * `"_WRM_type_ parameters`" ``[_WRM_type_-spec]`` See below for the required
     parameter spec for each type.

Example:

.. code-block:: xml

  <ParameterList name="PKs" type="ParameterList">
    ...
    <ParameterList name="flow" type="ParameterList">
      ...
      <ParameterList name="water retention evaluator" type="ParameterList">
        <Parameter name="minimum rel perm cutoff" type="double" value=" 0" />
        <Parameter name="use surface rel perm" type="bool" value="true" />
        <Parameter name="model parameters" type="string" value="WRM parameters" />
        ...
      </ParameterList>
      ...
    </ParameterList>
    ...
  </ParameterList>

  <ParameterList name="state" type="ParameterList">
    <ParameterList name="model parameters" type="ParameterList">
      <ParameterList name="WRM parameters" type="ParameterList">
        <ParameterList name="domain" type="ParameterList">
          <Parameter name="region" type="string" value="domain" />
          <Parameter name="wrm type" type="string" value="van Genuchten" />
          <Parameter name="van Genuchten alpha [Pa^-1]" type="double" value="2e-05" />
          <Parameter name="van Genuchten n [-]" type="double" value="1.58" />
          <Parameter name="residual saturation [-]" type="double" value="0.2" />
          <Parameter name="smoothing interval width [saturation]" type="double" value="0.05" />
          <Parameter name="dessicated zone thickness [m]" type="double" value="0.1" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
    ...
  </ParameterList>

*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_

#include "wrm_partition.hh"
#include "wrm.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
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
  static Utils::RegisteredFactory<Evaluator, WRMEvaluator> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
