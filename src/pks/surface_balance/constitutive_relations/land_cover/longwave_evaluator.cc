/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates incoming longwave radiation from vapor pressure and air temperature.

/*!

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec


    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"vapor pressure air key`" ``[string]`` **DOMAIN-vapor_pressure_air**

*/

#include "Key.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "longwave_evaluator.hh"


namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

LongwaveEvaluator::LongwaveEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  auto domain = Keys::getDomain(my_key_);
  air_temp_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(air_temp_key_);
  vp_air_key_ = Keys::readKey(plist, domain, "vapor pressure air", "vapor_pressure_air");
  dependencies_.insert(vp_air_key_);

  scale_ = plist.get<double>("scaling factor [-]", 1.0);
}

// Required methods from SecondaryVariableFieldEvaluator
void
LongwaveEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& air_temp = *S->GetFieldData(air_temp_key_)->ViewComponent("cell", false);
  const auto& vp_air = *S->GetFieldData(vp_air_key_)->ViewComponent("cell", false);
  auto& res = *result->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    res[0][c] = scale_ * Relations::IncomingLongwaveRadiation(air_temp[0][c], vp_air[0][c]);
  }
}

} //namespace
} //namespace
} //namespace

