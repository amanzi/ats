/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates relative permeability using water retention models.
/*!

Uses a list of regions and water retention models on those regions to evaluate
relative permeability, typically as a function of liquid saturation.

Most of the parameters are provided to the WRM model, and not the evaluator.
Typically these share lists to ensure the same water retention curves, and this
one is updated with the parameters of the WRM evaluator.  This is handled by
flow PKs.

Some additional parameters are available.

`"evaluator type`" = `"WRM rel perm`"

.. _rel-perm-evaluator-spec:
.. admonition:: rel-perm-evaluator-spec

   * `"use density on viscosity in rel perm`" ``[bool]`` **true**

   * `"boundary rel perm strategy`" ``[string]`` **boundary pressure** Controls
     how the rel perm is calculated on boundary faces.  Note, this may be
     overwritten by upwinding later!  One of:

      - `"boundary pressure`" Evaluates kr of pressure on the boundary face, upwinds normally.
      - `"interior pressure`" Evaluates kr of the pressure on the interior cell (bad idea).
      - `"harmonic mean`" Takes the harmonic mean of kr on the boundary face and kr on the interior cell.
      - `"arithmetic mean`" Takes the arithmetic mean of kr on the boundary face and kr on the interior cell.
      - `"one`" Sets the boundary kr to 1.
      - `"surface rel perm`" Looks for a field on the surface mesh and uses that.

   * `"minimum rel perm cutoff`" ``[double]`` **0.** Provides a lower bound on rel perm.

   * `"permeability rescaling`" ``[double]`` Typically rho * kr / mu is very big
     and K_sat is very small.  To avoid roundoff propagation issues, rescaling
     this quantity by offsetting and equal values is encourage.  Typically 10^7 or so is good.

   * `"WRM parameters`" ``[wrm-typedinline-spec-list]``  List (by region) of WRM specs.

   KEYS:

   - `"rel perm`"
   - `"saturation_liquid`"
   - `"density`" (if `"use density on viscosity in rel perm`" == true)
   - `"viscosity`" (if `"use density on viscosity in rel perm`" == true)
   - `"surface relative permeability`" (if `"boundary rel perm strategy`" == `"surface rel perm`")

*/

#ifndef AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

enum class BoundaryRelPerm {
  BOUNDARY_PRESSURE,
  INTERIOR_PRESSURE,
  HARMONIC_MEAN,
  ARITHMETIC_MEAN,
  ONE,
  SURF_REL_PERM
};

class RelPermEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit RelPermEvaluator(Teuchos::ParameterList& plist);
  RelPermEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  RelPermEvaluator(const RelPermEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrms_;
  Key sat_key_;
  Key dens_key_;
  Key visc_key_;
  Key surf_rel_perm_key_;

  bool is_dens_visc_;
  Key surf_domain_;
  BoundaryRelPerm boundary_krel_;

  double perm_scale_;
  double min_val_;

 private:
  static Utils::RegisteredFactory<Evaluator, RelPermEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
