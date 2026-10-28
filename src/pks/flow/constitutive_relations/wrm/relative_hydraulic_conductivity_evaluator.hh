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

  // Overridden (rather than relying on the generic per-dependency chain
  // rule) because, when use_surface_relperm_ is true, the "boundary_face"
  // component's contribution from krel_key_ must be replaced by
  // surf_krel_key_'s own derivative wrt the surface-domain counterpart of
  // wrt_key -- a cross-mesh substitution the generic chain rule cannot
  // perform, since it only ever asks dependencies about the SAME wrt_key
  // it was itself asked about. dens_key_ and visc_key_, and krel_key_'s
  // "cell" contribution, still go through the ordinary, generic chain rule.
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override;

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Requires surf_krel_key_'s own derivatives wrt the surface-domain
  // counterpart of every wrt_key this evaluator's own derivative has been
  // required for -- needed because the generic EnsureCompatibility_DepDerivs_
  // only requires a dependency's derivative wrt the SAME wrt_key, which is
  // never true for surf_krel_key_ (it depends on a surface-domain key).
  void EnsureCompatibility_SurfDepDerivs_(State& S);

 protected:
  static const std::string eval_type;

  bool use_surface_relperm_;
  KeyTag surf_krel_key_;
  Key surf_domain_; // only resolved/used when use_surface_relperm_

  KeyTag dens_key_, visc_key_, krel_key_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, RelativeHydraulicConductivityEvaluator> reg_;
};


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
