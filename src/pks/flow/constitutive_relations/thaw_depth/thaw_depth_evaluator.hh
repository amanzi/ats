/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An evaluator for calculating the depth to frozen soil/permafrost, relative to the surface.

/*!

Note that postive values are below the surface, and negative values
are not possible. 

Evaluator name: `"thaw depth`"

.. _thaw-depth-spec:
.. admonition:: thaw-depth-spec

    * `"transition width [K]`" ``[double]`` Width of the freeze curtain transition.

    KEYS:
      `"temperature`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ThawDepthEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit ThawDepthEvaluator(Teuchos::ParameterList& plist);
  ThawDepthEvaluator(const ThawDepthEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new ThawDepthEvaluator(*this));
  }

  // Custom HasFieldChanged forces this to be updated once.
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) override;
  // Custom EnsureCompatibility deals with multiple meshes
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;
  
 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override {}

 protected:    
  bool updated_once_;
  double trans_width_;
  Key domain_, domain_ss_;
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,ThawDepthEvaluator> reg_;

};
  
} //namespace
} //namespace 
