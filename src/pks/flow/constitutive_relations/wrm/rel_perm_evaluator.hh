/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Uses a list of regions and water retention models on those regions to evaluate
relative permeability, typically as a function of liquid saturation.

Most of the parameters are provided to the WRM model, and not the evaluator.
Typically these share lists to ensure the same water retention curves, and this
one is updated with the parameters of the WRM evaluator.  This is handled by
flow PKs.

Some additional parameters are available.

`"evaluator type`" = `"relative permeability, water retention model`"

.. _evaluator-relative-permeability-water-retention-model-spec:
.. admonition:: evaluator-relative-permeability-water-retention-model-spec

   * `"use density on viscosity in rel perm`" ``[bool]`` **true**

   * `"boundary rel perm strategy`" ``[string]`` **boundary pressure** Controls
     how the rel perm is calculated on boundary faces.  Note, this may be
     overwritten by upwinding later!  One of:

      - `"boundary pressure`" Evaluates kr of pressure on the boundary face, upwinds normally.
      - `"interior pressure`" Evaluates kr of the pressure on the interior cell (bad idea).
      - `"harmonic mean`" Takes the harmonic mean of kr on the boundary face and kr on the interior cell.
      - `"arithmetic mean`" Takes the arithmetic mean of kr on the boundary face and kr on the interior cell.
      - `"one`" Sets the boundary kr to 1.
      - `"surface rel perm`" Looks for a field on the surface mesh and uses that.

   * `"minimum rel perm cutoff`" ``[double]`` **0.** Provides a lower bound on rel perm.

   * `"model parameters`" ``[string]`` **"WRM parameters"** This will copy the
     named set of parameters, which must be in :ref:`State`'s `"model
     parameters`" list, into this list to evaluate relative
     permeability. Usually this shares a common set of parameters (usually
     called `"WRM parameters`"), with the :ref:`WRM Evaluator`.  Some models
     may use a customed name, e.g., `"relative permeability parameters`", and
     declare `"relative permeability parameters`" in `"model parameters`" under
     state, which allows one to use different models for the WRM and relative
     permeability.

   KEYS:

   - `"relative permeability`"

   DEPENDENCIES:

   - `"saturation_liquid`"
   - `"density`" (if `"use density on viscosity in rel perm`" == true)
   - `"viscosity`" (if `"use density on viscosity in rel perm`" == true)
   - `"surface relative permeability`" (if `"boundary rel perm strategy`" == `"surface rel perm`")



**Example 1:**

Using the same set of van Genuchten model paramters for WRM and relative permeability

.. code-block:: xml

  <ParameterList name="PKs" type="ParameterList">
    ...
    <ParameterList name="flow" type="ParameterList">
      ...
      <ParameterList name="water retention evaluator" type="ParameterList">
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
    <ParameterList name="evaluators" type="ParameterList">
      <ParameterList name="relative_permeability" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="relative permeability, van Genuchten" />
        <Parameter name="model parameters" type="string" value="WRM parameters" />
        <Parameter name="use surface rel perm" type="bool" value="true" />
        <Parameter name="minimum rel perm cutoff" type="double" value=" 0" />
      </ParameterList>
      ...
    </ParameterList>
    ...
  </ParameterList>


**Example 2:**

Using different set of model/paramters for WRM and relative permeability,
van Genuchten for WRM and Brooks-Corey for relative permeability. Note that
in this case, van Genuchten parameters and Brooks-Corey parameters need to
be consistent. Using tool `"convert_parameters_vg2bc.py`" to convert van Genuchten
parameters to Brooks-Corey parameters.

.. code-block:: xml

  <ParameterList name="PKs" type="ParameterList">
    ...
    <ParameterList name="flow" type="ParameterList">
      ...
      <ParameterList name="water retention evaluator" type="ParameterList">
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
      <ParameterList name="relative permeability parameters" type="ParameterList">
        <ParameterList name="domain" type="ParameterList">
          <Parameter name="region" type="string" value="domain" />
          <Parameter name="wrm type" type="string" value="Brooks-Corey" />
          <Parameter name="Brooks-Corey lambda [-]" type="double" value="0.49" />
          <Parameter name="Brooks-Corey saturated matric suction [Pa]" type="double" value="32439.03" />
          <Parameter name="residual saturation [-]" type="double" value="0.2" />
          <Parameter name="smoothing interval width [saturation]" type="double" value="0.05" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
    <ParameterList name="evaluators" type="ParameterList">
      <ParameterList name="relative_permeability" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="relative permeability, van Genuchten" />
        <Parameter name="model parameters" type="string" value="relative permeability parameters" />
        <Parameter name="use surface rel perm" type="bool" value="true" />
        <Parameter name="minimum rel perm cutoff" type="double" value=" 0" />
      </ParameterList>
      ...
    </ParameterList>
    ...
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {      
namespace Flow {

enum class BoundaryRelPerm {
  BOUNDARY_PRESSURE,
  INTERIOR_PRESSURE,
  HARMONIC_MEAN,
  ARITHMETIC_MEAN,
  ONE,
  SURF_REL_PERM
};

class RelPermEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit RelPermEvaluator(Teuchos::ParameterList& plist);
  RelPermEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  RelPermEvaluator(const RelPermEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrms_;
  Key sat_key_;
  Key dens_key_;
  Key visc_key_;
  Key surf_rel_perm_key_;

  bool is_dens_visc_;
  Key surf_domain_;
  BoundaryRelPerm boundary_krel_;

  double min_val_;

 private:
  static Utils::RegisteredFactory<Evaluator, RelPermEvaluator> factory_;
};

} // namespace Flow
}  
} // namespace Amanzi

#endif
