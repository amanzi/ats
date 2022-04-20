/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Plant wilting factor provides a moisture availability-based limiter on transpiration.
/*!

  Generated via evaluator_generator with:
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

Note this makes use of LandCover objects for mafic potential of fully open and
fully closed stomata.

Note the challenges of using this model with arbitrary van Genuchten WRMs.  See
Verhoef & Egea, Ag. & Forest Meteorology, 2014
https://doi.org/10.1016/j.agrformet.2014.02.009


.. _plant-wilting-factor-evaluator-spec:
.. admonition:: plant-wilting-factor-evaluator-spec

   KEYS:
   - `"capillary pressure`" **DOMAIN-capillary_pressure_gas_liq**


*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class PlantWiltingFactorModel;

class PlantWiltingFactorEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PlantWiltingFactorEvaluator(Teuchos::ParameterList& plist);
  PlantWiltingFactorEvaluator(const PlantWiltingFactorEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key pc_key_;

  Key domain_surf_;
  Key domain_sub_;

  LandCoverMap land_cover_;
  std::map<std::string,Teuchos::RCP<PlantWiltingFactorModel>> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PlantWiltingFactorEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

