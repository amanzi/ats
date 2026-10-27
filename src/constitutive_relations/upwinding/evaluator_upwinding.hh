/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Base class for evaluators that upwind a cell-based field onto faces.
/*!

Upwinding takes a cell-based coefficient (typically a relative permeability or
relative conductivity) and produces a face-based value, used by a div-grad
operator's nonlinear coefficient.  Rules differ in how they decide, per face,
which neighboring cell's value (or what blend of the two) to use -- but all
share the same output structure (FACE, 1 DoF for the value) and the same
derivative structure (FACE, 2 DoF for the derivative).

The 2 DoFs are per PRESSURE, not per cell: slot 0 is always d(value)/d(p0)
and slot 1 is always d(value)/d(p1), where {p0, p1} are the two pressures the
face value can depend on. For an INTERIOR face, both are neighboring cells'
own pressures, ordered by the mesh's own `getFaceCells(f)` convention. For a
BOUNDARY face, p0 is the one real neighbor cell's pressure and p1 is the
boundary candidate's own pressure. Both slots are always computed as real,
separate numbers -- never collapsed into one slot by upwind direction, even
in the hard-upwind case (where the non-governing side's own slot is a
genuine, computed zero, not an absent/repurposed slot).

The boundary candidate is always `field_key_`'s own `"boundary_face"`
component -- this evaluator has no knowledge of where that value came from
(e.g. the ordinary WRM-at-boundary-pressure chain, or a surface-domain
override stamped in by another evaluator upstream, such as
RelativeHydraulicConductivityEvaluator's `"use surface rel perm"`). It only
distinguishes two upwind POLICIES, selected by `"clobber policy"`:

- `"none"` (default): blend the interior cell's value against
  `field_key_`'s own `"boundary_face"` value using the ordinary per-rule
  upwind decision, exactly as at an interior face (just with one "cell"
  replaced by the boundary_face component).
- `"clobber"`: skip the upwind decision entirely -- always use
  `field_key_`'s own `"boundary_face"` value, regardless of flux direction.

Because the dependencies of this evaluator live on different mesh entity
kinds than its own output (a CELL-based field to upwind, plus a directional
dependency that may be CELL or FACE depending on the rule), this evaluator
cannot use EvaluatorSecondaryMonotypeCV's default same-structure
EnsureCompatibility hooks.  Those are overridden here; the generic/type-only
hooks (claim ownership, flags, deps existence, claim/dep derivs, recursing
into dependency EnsureCompatibility, units) are inherited unmodified.

KEYS:

- `"field key`" The cell-based field to be upwinded.

*/

#ifndef AMANZI_UPWINDING_EVALUATOR_UPWINDING_HH_
#define AMANZI_UPWINDING_EVALUATOR_UPWINDING_HH_

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class EvaluatorUpwinding : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaluatorUpwinding(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorUpwinding(const EvaluatorUpwinding& other) = default;

 protected:
  // Structure: my value is on FACE, 1 DoF.
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // Structure: my derivative is on FACE, 2 DoF -- NOT copied from my own
  // value structure, since the value is 1 DoF and the derivative is 2 DoF.
  virtual void EnsureCompatibility_DerivStructure_(State& S) override;

  // Structure: push explicit, named, per-dependency requirements instead of
  // the CV default (which pushes a single factory uniformly to every
  // dependency).  The field-to-upwind dependency always needs CELL plus
  // BOUNDARY_FACE; each concrete rule adds its own requirement for its
  // directional dependency in its own override, calling this one first.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Fully custom chain rule: pulls d(field)/dp on cells, then attributes it
  // to the 2-DoF face-based derivative using the same per-rule upwind
  // decision that Evaluate_ uses for the value.  Does NOT use the inherited
  // elementWiseMultiply chain-rule loop, which assumes matching structure
  // between the derivative and its dependencies' derivatives -- false here.
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override;

  // Per-rule upwind decision + value, called from Update_ via the base
  // class's Evaluate_ dispatch.
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override = 0;

  // Only ever has one real branch (the field-to-upwind dependency) -- the
  // directional dependency is never differentiated through, so this is
  // never asked to differentiate wrt it (see IsDifferentiableWRT filtering
  // in the base EvaluatorSecondaryMonotype<CV>::UpdateDerivative_ loop).
  // Not used by this class's own custom UpdateDerivative_ (which calls
  // UpwindDerivative_ instead), but kept pure virtual to satisfy the base
  // class's abstract interface for any caller going through the generic
  // path.
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override = 0;

  // Per-rule: given cell-based d(field)/dp, fill the 2-DoF face-based
  // derivative CV using the same upwind decision as Evaluate_.
  virtual void UpwindDerivative_(const State& S,
                                 const CompositeVector& dfield_dp,
                                 CompositeVector& result) const = 0;

  // Seeds boundary/exterior faces of `result` from field_key_'s own
  // "boundary_face" component.  This must be called BEFORE the per-rule
  // upwind loop runs: a boundary face has no "other side" cell, so the
  // upwind decision for such faces (in particular the near-zero-flux
  // blended case) reads back whatever is already sitting in `result`'s
  // face component as the missing side's value.  Without this seed, that
  // read is of uninitialized/stale data. A no-op when field_key_ has no
  // "boundary_face" component at all.
  //
  // When is_derivative is true, `wrt_key`/`wrt_tag` select field_key_'s own
  // "boundary_face" derivative wrt wrt_key to seed instead of its value.
  void SeedBoundaryFaces_(const State& S,
                          CompositeVector& result,
                          bool is_derivative,
                          const Key& wrt_key = "",
                          const Tag& wrt_tag = Tags::DEFAULT) const;

  // "Clobber" policy.  After the per-rule Evaluate_/UpwindDerivative_ core
  // logic runs, if strategy_ == "clobber", unconditionally OVERWRITE
  // boundary faces of `result` again with field_key_'s own "boundary_face"
  // value, undoing whatever the upwind decision chose there. This reuses
  // the same seeding logic as SeedBoundaryFaces_, just called a second
  // time, after the upwind loop.
  void ApplyClobber_(const State& S,
                     CompositeVector& result,
                     bool is_derivative,
                     const Key& wrt_key = "",
                     const Tag& wrt_tag = Tags::DEFAULT) const;

 protected:
  Key field_key_;
  std::string strategy_; // "none" or "clobber"
};

} // namespace Relations
} // namespace Amanzi

#endif
