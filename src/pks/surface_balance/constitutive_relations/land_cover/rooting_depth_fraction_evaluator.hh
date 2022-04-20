/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Provides a depth-based profile of root density.
/*!

Sets the root fraction as a function of depth,

.. math:
   F_root =  ( \alpha \; exp(-\alpha z) + \beta \; exp(-\beta z) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

Note that all three parameters, a, b, and the cutoff, are provided in the
LandCover type.

.. _rooting-depth-fraction-evaluator-spec:
.. admonition:: rooting-depth-fraction-evaluator-spec

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Sane default provided for most domain names.

   KEYS:
   - `"depth`" **DOMAIN-depth**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **SURFACE_DOMAIN-cell_volume**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RootingDepthFractionModel;

class RootingDepthFractionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  RootingDepthFractionEvaluator(Teuchos::ParameterList& plist);
  RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  void InitializeFromPlist_();

  Key z_key_;
  Key cv_key_;
  Key surf_cv_key_;

  Key domain_surf_;
  Key domain_sub_;

  LandCoverMap land_cover_;
  std::map<std::string, Teuchos::RCP<RootingDepthFractionModel>> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RootingDepthFractionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

