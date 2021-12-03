/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Ahmad Jan
*/
//! Calculates the change in surface elevation of a column mesh over time.

/*!
.. _subsidence-evaluator-spec
.. admonition:: subsidence-evaluator-spec

   DEPENDENCIES:
   - `"initial elevation`" **DOMAIN-initial_elevation** Elevation at the start of the simulation.
   - `"base porosity`" **DOMAIN_SUBSURF-base_porosity** A target that changes
     when the mesh changes -- this is a workaround for the fact that the mesh
     change is not in the DAG.

NOTE: this only works on columns and should go away!  It may be the case that
this can be replaced by the typical overland flow elevation, but that may not
be true if it is not consistent with the face_centroid value of the mesh.
Some thought may be required before making that substitution.
*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class SubsidenceEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SubsidenceEvaluator(Teuchos::ParameterList& plist);
  SubsidenceEvaluator(const SubsidenceEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override {}

  bool updated_once_;
  Key domain_;
  Key bp_key_, init_elev_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubsidenceEvaluator> reg_;

};

} //namespace
} //namespace

