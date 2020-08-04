/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (coonet @ ornl.gov)
   
 ------------------------------------------------------------------------- */

//! AlbedoSubgridEvaluator: evaluates albedos and emissivities with a subgrid model.

/*!
Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for three channels -- water/ice, land, and
snow.  Note this internally calculates albedo of snow based upon snow density.

Channels are: 0 = land, 1 = water/ice, 2 = snow.

* `"albedo ice [-]`" ``[double]`` **0.44** 
* `"albedo water [-]`" ``[double]`` **0.1168** 
* `"albedo ground surface [-]`" ``[double]`` **0.135** Defaults to that of tundra.

* `"emissivity ice [-]`" ``[double]`` **0.98** 
* `"emissivity water [-]`" ``[double]`` **0.995** 
* `"emissivity ground surface [-]`" ``[double]`` **0.92** Defaults to that of tundra.
* `"emissivity snow [-]`" ``[double]`` **0.98**

* `"snow density key`" ``[string]`` **SNOW_DOMAIN-density** 
* `"unfrozen fraction key`" ``[string]`` **DOMAIN-unfrozen_fraction**

*/


#ifndef ALBEDO_SUBGRID_EVALUATOR_HH_
#define ALBEDO_SUBGRID_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class AlbedoSubgridEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  AlbedoSubgridEvaluator(Teuchos::ParameterList& plist);
  AlbedoSubgridEvaluator(const AlbedoSubgridEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AlbedoSubgridEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_;
  Key unfrozen_fraction_key_;

  double a_water_, a_ice_, a_tundra_;
  double e_water_, e_ice_, e_tundra_, e_snow_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,AlbedoSubgridEvaluator> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
