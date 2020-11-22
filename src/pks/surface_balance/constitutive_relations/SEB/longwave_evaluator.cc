/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates incoming longwave radiation from rel humidity and air temperature.

/*!

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec

    * `"minimum relative humidity [-]`" ``[double]`` **0.1** Sets a minimum rel humidity, RH=0 breaks the model.

    DEPENDENCIES:
    
    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"relative humidity key`" ``[string]`` **DOMAIN-relative_humidity**
         
*/

#include "Key.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "longwave_evaluator.hh"


namespace Amanzi {
namespace SurfaceBalance {

LongwaveEvaluator::LongwaveEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  auto domain = Keys::getDomain(my_key_);
  air_temp_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(air_temp_key_);
  rel_hum_key_ = Keys::readKey(plist, domain, "relative humidity", "relative_humidity");
  dependencies_.insert(rel_hum_key_);

  min_rel_hum_ = plist.get<double>("minimum relative humidity [-]", 0.1);  
  stephB_ = SEBPhysics::ModelParams().stephB;
}

// Required methods from SecondaryVariableFieldEvaluator
void
LongwaveEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& air_temp = *S->GetFieldData(air_temp_key_)->ViewComponent("cell", false);
  const auto& rel_hum = *S->GetFieldData(rel_hum_key_)->ViewComponent("cell", false);
  auto& res = *result->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    res[0][c] = SEBPhysics::CalcIncomingLongwave(air_temp[0][c], std::max(min_rel_hum_, rel_hum[0][c]), stephB_);
  }  
}

} //namespace
} //namespace

