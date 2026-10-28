# TODO: Epetra -> Tpetra port

Tracking list for structural gaps in the Epetra -> Tpetra port, as opposed to
physics-specific porting work (individual PKs, evaluators, etc).

## 0. Update all coupled-water tests with a non-injected MPC analogue, compare solutions

For every existing "coupled water" (old, injected `MPCCoupledWater`) test, add
a `MPCCoupledWaterNonInjected` analogue and compare the actual solutions
between the two (not just that both run) -- similar to what was already done
for the `column_drainage`/`column_exfiltration`/etc. family this session.
Planned to be picked up ASAP (next session or so).

## 1. Overland flow: several Manning upwinding schemes were dropped, not ported

The Epetra-era `src/operators/upwinding/` directory (on `ecoon/io`) implemented
a factory (`UpwindFluxFactory`, selected via overland flow's `"upwind type"`
plist parameter) with six named schemes. Of these, only `"manning upwind"`
(-> `EvaluatorUpwindTotalFlux`) and `"manning cell centered"` (a local
fallback) are supported by the current `OverlandPressureFlow` PK. The
following were genuinely dropped during the port, not just left unported in a
working state:

- `"manning harmonic mean"` (`UpwindFluxHarmonicMean`) -- NOTE: the old code's
  own factory comment flagged this as dangerous/deprecated ("this option has
  shown up in a few input files and broken user code") but left the option
  live anyway; worth deciding intentionally whether to port it at all.
- `"manning split denominator"` (`UpwindFluxSplitDenominator`) -- a genuine,
  still-useful Manning-flow scheme using slope/Manning-coefficient fields.
- `"manning ponded depth passthrough"` (`UpwindFluxFOCont`) -- same, another
  genuine Manning-flow scheme.

(`"manning elevation stabilized"` / `UpwindElevationStabilized` was already
fully commented out of the old factory on `ecoon/io` -- dead before the port,
nothing lost there.)

`UpwindFluxFactory` itself does not exist in the ported code at all; the stale
file `src/pks/flow/overland_pk.cc` (excluded from the build via
`src/pks/flow/CMakeLists.txt`) still `#include`s the deleted
`UpwindFluxFactory.hh` -- harmless while that file stays out of the build, but
should be cleaned up (ported or deleted) rather than left as silent dead code.

No current XML config requests anything other than the default
`"manning upwind"`, so nothing is live-broken today -- but any new overland
flow config wanting one of the dropped schemes will hit a hard
"not supported in this build" error.

Separately: `EvaluatorUpwindPotentialDifference` (a *different* scheme, used
by Richards/subsurface flow rather than overland) WAS ported, compiles,
links, and is registered under `"upwind: potential difference"` -- just not
exercised by any current test. No action needed there.

## 2. `EnsureEvaluators` is not called early enough to make `IsDifferentiableWRT`/`HasDependency`/`ProvidesKey` safe during PK `setup()`

On `master`, `State::RequireEvaluator(key, tag)` is the single chokepoint for
obtaining/constructing an evaluator, and it calls the newly-constructed
evaluator's own `EnsureEvaluators(*this)` immediately, inline (see
`~/code/ats/ats/repos/master/src/state/State.cc`, `RequireEvaluator`, calls at
the "provided by another evaluator" branch and the "create from plist"
branch). `EnsureEvaluators` recursively `RequireEvaluator`s every dependency,
which in turn (for any dependency that is itself newly constructed) triggers
its own `EnsureEvaluators` the same way. The net effect: the instant any code
calls `RequireEvaluator` on a key, the ENTIRE dependency subtree rooted at
that evaluator is guaranteed fully constructed and valid -- not deferred.
This means a PK's own `setup()` can safely call `IsDifferentiableWRT`,
`HasDependency`, or `ProvidesKey` on any evaluator right after requiring it,
with no further caveats.

On this (tpetra) branch, `EnsureEvaluators` exists on `State` and on the
`Evaluator` base classes (same virtual method, same recursive-dependency
logic) but is only invoked from ONE place: a single global pass inside
`State::Setup()` (`src/state/State.cc`), which walks the (by-then-complete)
`evaluators_` map once, AFTER every PK in the PK tree has already run its own
`setup()`. `RequireEvaluator` itself does NOT call `EnsureEvaluators` on this
branch. Consequence: during any individual PK's `setup()`/`parseParameterList()`
-- which is exactly when a PK would want to check `IsDifferentiableWRT` before
deciding whether to `RequireDerivative` on some coefficient -- the dependency
subtree for a freshly-required evaluator is NOT yet guaranteed valid. Calling
`IsDifferentiableWRT` at that point can give a wrong (or even
crashing/assert-failing) answer depending on construction order.

This is why `Richards::SetupRichardsFlow_`'s `IsDifferentiableWRT` guard
around `RequireDerivative(coef_key_, ...)` had to stay commented out/removed
this session (with the original code's own comment: "this is not valid until
State always EnsureEvaluators! --ETC") -- the one place it's actually SAFE to
call `IsDifferentiableWRT` on this branch is deep in the timestepping loop
(`UpdatePermeabilityDerivativeData_`), well after `State::Setup()` has
completed, never inside `setup()` itself.

**Fix**: port master's `RequireEvaluator`-calls-`EnsureEvaluators`-inline
design to this branch's `State::RequireEvaluator`. Likely requires auditing
PK `parseParameterList()`/`setup()` ordering too (per master's convention,
primary variables' evaluator plists get pushed into state's `"evaluators"`
list during `parseParameterList()`, before any secondary evaluator that
depends on them gets constructed, so the dependency leaves exist by the time
anything needs them) -- confirm this branch's PKs already follow that
convention, or fix them to, as part of this change.

## 3. Audit: systematic pass through Amanzi/ATS `master` history for other structural changes not yet ported

Item 2 above was found by noticing one specific symptom (a disabled
`IsDifferentiableWRT` guard with an explanatory comment) and tracing it back
to a real, intentional design change on `master` that this branch never
picked up. There are likely other structural/framework-level changes on
`master` from the last 1-2 years (state/evaluator framework, PK lifecycle
conventions, etc. -- as opposed to individual PK physics, which this branch
has been actively porting PK-by-PK) that would be worth finding and either
porting now or consciously deferring, before sinking more effort into
physics-level porting on top of a framework that may still be missing other
master-side improvements.

Not yet scoped (which directories, what date range, what counts as
"structural" vs. "physics") -- do that scoping when this item is actually
picked up, likely as its own dedicated session/pass rather than an aside
during physics work.

## 4. New idea: column-structured preconditioner/solver for device

Not a porting gap -- a new design, doesn't exist anywhere in `master` or this
branch today. Many ATS meshes (columnar extruded 2D-surface x 1D-vertical
meshes, the common case for integrated hydrology problems) have a natural
per-column structure: within a single column, cells/faces are contiguous and
ordered top-to-bottom, and the vertical (z-direction) coupling along a column
is tridiagonal (or near-tridiagonal) in a cell-centered or face-based
discretization. Idea: exploit this directly rather than treating the whole
mesh as a generic unstructured sparse system:

- Lay out each column's unknowns contiguously in memory (column-major
  ordering across the whole mesh, not just within a column).
- Solve the vertical (within-column) coupling exactly and cheaply via a
  tridiagonal solve (e.g. Thomas algorithm) -- O(n) per column, no fill-in,
  no generic sparse solver overhead.
- Columns are independent of each other for this vertical solve, so this
  parallelizes trivially and efficiently across columns on device (one
  column per thread/warp, contiguous memory access pattern).
- Horizontal (cross-column) coupling is the remaining piece -- handle via
  some outer iterative method (e.g. a column-block Jacobi/Gauss-Seidel, or
  using the per-column tridiagonal solve as a preconditioner for an outer
  Krylov method on the horizontal coupling only).

Motivation: likely a much better fit for device/Kokkos parallelism than a
generic AMG/ILU preconditioner over the whole unstructured mesh, for the
large fraction of ATS problems that are actually columnar. Not yet scoped
beyond this description -- no design doc, no prototype, no decision on where
this would live (a new `PDE_Diffusion` discretization? a standalone
`Operators::Preconditioner` implementation? something in `solvers/`?).
