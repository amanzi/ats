/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An evaluator for calculating the depth to water table, relative to the surface.

/*!

Note that postive values are below the surface, and negative values
are not possible.  One may choose to, in post-processing, do:

  water_table_height = where(ponded_depth > 0, ponded_depth, -water_table_depth)

to calculate a water table height that is valid for both surface and
subsurface water.

Evaluator name: `"water table depth`"

.. _water-table-depth-spec:
.. admonition:: water-table-depth-spec

    KEYS:
      `"saturation_gas`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class WaterTableDepthEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit WaterTableDepthEvaluator(Teuchos::ParameterList& plist);
  WaterTableDepthEvaluator(const WaterTableDepthEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new WaterTableDepthEvaluator(*this));
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
  Key sat_key_;
  Key domain_, domain_ss_;

 private:
  static Utils::RegisteredFactory<Evaluator,WaterTableDepthEvaluator> reg_;

};
  
} //namespace
} //namespace 
