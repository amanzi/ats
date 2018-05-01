/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An evaluator for converting the darcy flux to volumetric flux

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/
#ifndef AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_
#define AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class Volumetric_FluxEvaluator :
    public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  explicit
  Volumetric_FluxEvaluator(Teuchos::ParameterList& plist);

  Volumetric_FluxEvaluator(const Volumetric_FluxEvaluator& other);

  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility(State& S) override;

  Teuchos::RCP<Evaluator> Clone() const;

protected:
  // Required methods from SecondaryVariableEvaluator
  virtual void Evaluate_(const State& S,
                         CompositeVector& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key, const Key& wrt_tag, CompositeVector& result) override;

  Key flux_key_, flux_tag_key_;
  Key dens_key_, dens_tag_key_;
  Key mesh_key_, mesh_tag_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,Volumetric_FluxEvaluator> fac_;

};


}//namespace
}//namespace

#endif
