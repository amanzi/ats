/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! PresElevEvaluator: evaluates h + z

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
Evaluator type: "snow skin potential"

.. math::
  h + z + h_{{snow}} + dt * P_{{snow}}

* `"my key`" ``[string]`` **snow_skin_potential** Names the potential variable evaluated [m]
* `"ponded depth key`" ``[string]`` **ponded_depth** Names the surface water depth variable. [m]
* `"snow depth key`" ``[string]`` **snow_depth** Names the snow depth variable. [m]
* `"precipitation snow key`" ``[string]`` **precipitation_snow** Names the snow precipitation key. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]
* `"dt factor [s]`" ``[double]`` A free-parameter factor for providing a time scale for diffusion of snow precipitation into low-lying areas.  Typically on the order of 1e4-1e7. This timestep times the wave speed of snow provides an approximate length of how far snow precip can travel.  Extremely tunable! [s]

NOTE: This is equivalent to a generic AdditiveEvaluator_

Example:

.. code-block:: xml

  <ParameterList name="snow_skin_potential" type="ParameterList">
    <Parameter name="field evaluator type" type="string" value="snow skin potential" />
    <Parameter name="dt factor [s]" type="double" value="864000.0" />
  </ParameterList>

*/


#ifndef AMANZI_FLOWRELATIONS_SNOW_SKIN_POTENTIAL_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SNOW_SKIN_POTENTIAL_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class SnowSkinPotentialEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SnowSkinPotentialEvaluator(Teuchos::ParameterList& plist);

  SnowSkinPotentialEvaluator(const SnowSkinPotentialEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Key precip_key_;
  Key sd_key_;
  Key pd_key_;
  Key elev_key_;
  double factor_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SnowSkinPotentialEvaluator> factory_;

};

} //namespace
} //namespace

#endif
