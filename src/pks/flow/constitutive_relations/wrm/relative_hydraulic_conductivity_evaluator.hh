/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

This computes the product,

.. math::
   k = \frac{n}{\mu} k_r

which is the scalar portion of the hydraulic conductance (excludes the absolute
permeability).  Relative permeability, k_r, is computed using a water retention
model, which provides generic methods for this quantity based on internal
parameters (e.g. van Genuchten/Mualem, Brooks & Corey, etc.).

This evaluator layers on a few details relative to the above equation,
including the ability to import a surface relative permeability for incoming
infiltration.

type : `"relative hydraulic conductivity`"

.. _relative-hydraulic-conductivity-evaluator-spec
.. admonition:: relative-hydraulic-conductivity-evaluator-spec

   * `"use surface rel perm`" ``[bool]]`` **false**                
   * `"use surface rel perm`" ``[bool]]`` **false**                

   KEYS
   - `"relative_permeability`" **DOMAIN-relative_permeability**
   - `"density`" **DOMAIN-molar_density_liquid**
   - `"viscosity`" **DOMAIN-viscosity**

*/

#pragma once

#include "Key.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class RelativeHydraulicConductivityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  RelativeHydraulicConductivityEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  RelativeHydraulicConductivityEvaluator(const RelativeHydraulicConductivityEvaluator& other) =
    default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  std::string getType() const override { return eval_type; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  static const std::string eval_type;

  bool use_surface_relperm_;
  KeyTag surf_krel_key_;

  KeyTag dens_key_, visc_key_, krel_key_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, RelativeHydraulicConductivityEvaluator> reg_;
};


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
