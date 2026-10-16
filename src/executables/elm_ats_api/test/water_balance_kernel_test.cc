/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rich Fiorella, Ethan Coon
*/

//! Unit tests for the ELM water-balance kernel used by ATS.
/*!

  Exercises the dependency-free, C-bound kernel elm_ats_water_balance_error_c
  (implemented in WaterBalanceKernel.F90 and declared in
  elm_balance_interface_private.hh) that ATS calls at the coupling boundary to
  gate step acceptance (see ELM_ATSDriver::checkELMWaterBalance_).

  The kernel computes ELM's column water-balance error [mm]:

    errh2o = (endwb - begwb)
             - (source - evap - tran - runoff - baseflow) * dtime

  where the C-facing subset maps to ELM's full formula as:
    source -> forc_rain (+), evap+tran -> qflx_evap_tot (sink),
    runoff -> qflx_surf (sink), baseflow -> qflx_drain (sink),
    dStorage = endwb - begwb.  All ELM-only terms are zero.

  Units: endwb/begwb [mm], all fluxes [mm s-1], dtime [s], errh2o [mm].

*/

#include <cmath>
#include <cstdio>
#include <vector>

#include "elm_balance_interface_private.hh"

namespace {

// --- Minimal dependency-free test harness ----------------------------------

int g_checks = 0;
int g_failures = 0;

// Records a floating-point comparison; prints a diagnostic and flags failure
// when |expected - actual| exceeds tol.
void
check_close(double expected, double actual, double tol, const char* file,
            int line, const char* expr)
{
  ++g_checks;
  const double diff = std::fabs(expected - actual);
  if (diff > tol) {
    ++g_failures;
    std::fprintf(stderr,
                 "%s:%d: FAILED CHECK_CLOSE(%s)\n"
                 "    expected %.17g, got %.17g (|diff|=%.3g > tol=%.3g)\n",
                 file, line, expr, expected, actual, diff, tol);
  }
}

#define CHECK_CLOSE(expected, actual, tol)                                     \
  check_close((expected), (actual), (tol), __FILE__, __LINE__, #actual)

// Independent C++ evaluation of the same formula, used as the reference oracle.
double
reference(double endwb, double begwb, double source, double evap, double tran,
          double baseflow, double runoff, double dtime)
{
  const double net_flux = source - evap - tran - runoff - baseflow;
  return (endwb - begwb) - net_flux * dtime;
}

constexpr double TIGHT = 1.0e-12;

// --- Test cases -------------------------------------------------------------

// Storage change exactly matches net flux over the step -> zero error.
void
PerfectBalance()
{
  const double dt = 1800.0;   // 30 min
  const double source = 2.0e-3;   // mm/s in
  const double begwb = 100.0;     // mm
  const double endwb = begwb + source * dt; // storage rises by the source input
  CHECK_CLOSE(0.0, elmWaterBalanceError(endwb, begwb, source, 0.0, 0.0, 0.0, 0.0, dt), TIGHT);
}

// Pure source, no storage change: all the input shows up as imbalance.
void
PureSource()
{
  const double dt = 1800.0;
  const double source = 2.0e-3;
  CHECK_CLOSE(-source * dt, elmWaterBalanceError(50.0, 50.0, source, 0.0, 0.0, 0.0, 0.0, dt), TIGHT);
}

// Each sink term, in isolation, with no storage change, adds +term*dt to error.
void
EvapSink()
{
  const double dt = 1800.0, evap = 1.5e-3;
  CHECK_CLOSE(evap * dt, elmWaterBalanceError(50.0, 50.0, 0.0, evap, 0.0, 0.0, 0.0, dt), TIGHT);
}

void
TranSink()
{
  const double dt = 1800.0, tran = 1.1e-3;
  CHECK_CLOSE(tran * dt, elmWaterBalanceError(50.0, 50.0, 0.0, 0.0, tran, 0.0, 0.0, dt), TIGHT);
}

void
RunoffSink()
{
  const double dt = 1800.0, runoff = 0.7e-3;
  CHECK_CLOSE(runoff * dt, elmWaterBalanceError(50.0, 50.0, 0.0, 0.0, 0.0, 0.0, runoff, dt), TIGHT);
}

void
BaseflowSink()
{
  const double dt = 1800.0, baseflow = 0.4e-3;
  CHECK_CLOSE(baseflow * dt, elmWaterBalanceError(50.0, 50.0, 0.0, 0.0, 0.0, baseflow, 0.0, dt), TIGHT);
}

// evap and tran both map to qflx_evap_tot, i.e. they are summed.
void
EvapPlusTran()
{
  const double dt = 1800.0, evap = 1.5e-3, tran = 1.1e-3;
  CHECK_CLOSE((evap + tran) * dt,
              elmWaterBalanceError(50.0, 50.0, 0.0, evap, tran, 0.0, 0.0, dt), TIGHT);
}

// dtime == 0: the flux term drops out; the kernel has no dt guard (that lives
// in the driver), so the result is exactly the storage change.
void
StorageOnlyZeroDt()
{
  CHECK_CLOSE(12.5 - 10.0,
              elmWaterBalanceError(12.5, 10.0, 3.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 0.0),
              TIGHT);
}

// All terms nonzero: kernel must match an independent evaluation of the formula.
void
CombinedNonZero()
{
  const double dt = 900.0;
  const double endwb = 123.456, begwb = 120.0;
  const double source = 3.3e-3, evap = 1.2e-3, tran = 0.9e-3;
  const double baseflow = 0.5e-3, runoff = 0.6e-3;
  CHECK_CLOSE(reference(endwb, begwb, source, evap, tran, baseflow, runoff, dt),
              elmWaterBalanceError(endwb, begwb, source, evap, tran, baseflow, runoff, dt),
              TIGHT);
}

// --- Unit conversions: ATS [mol] storage -> ELM [mm] column storage ----------
//
// These cover the arithmetic ELM_ATSDriver::checkELMWaterBalance_ performs
// before calling the kernel, which the cases above do not touch because they
// feed mm values in directly.  The driver passes Epetra row pointers and a mesh
// column view; here plain arrays and a vector stand in for them.

// Single column, non-uniform molar density: volume is the sum of per-cell
// mol/(mol m-3), so a hoisted or swapped divisor changes the answer.
void
ColumnVolumeNonUniformDensity()
{
  // 3 subsurface cells + the surface cell, all with distinct molar densities.
  const double wc[3] = { 5.5e4, 4.0e4, 3.0e4 };       // mol
  const double n_liq[3] = { 5.5e4, 5.0e4, 6.0e4 };    // mol m-3
  const std::vector<int> cells = { 0, 1, 2 };
  const double surf_wc = 1.0e3, surf_n_liq = 5.0e4;   // mol, mol m-3

  // 1.0e3/5.0e4 + 5.5e4/5.5e4 + 4.0e4/5.0e4 + 3.0e4/6.0e4
  //   = 0.02 + 1.0 + 0.8 + 0.5 = 2.32 m^3
  CHECK_CLOSE(2.32, elmColumnWaterVolume(surf_wc, surf_n_liq, wc, n_liq, cells), TIGHT);
}

// The column's cell ids need not be contiguous or ordered: ELM columns index
// into a full subsurface array.  Only the listed cells may contribute -- a
// version that walked 0..n-1 instead of the id list would fail here.
void
ColumnVolumeScatteredCellIds()
{
  // 8-cell subsurface array; this column owns only {7, 2, 5}.
  const double wc[8] = { 9.9e9, 9.9e9, 2.0e4, 9.9e9, 9.9e9, 3.0e4, 9.9e9, 1.0e4 };
  const double n_liq[8] = { 1.0, 1.0, 4.0e4, 1.0, 1.0, 5.0e4, 1.0, 2.0e4 };
  const std::vector<int> cells = { 7, 2, 5 };
  const double surf_wc = 0.0, surf_n_liq = 5.0e4;

  // 0 + 1.0e4/2.0e4 + 2.0e4/4.0e4 + 3.0e4/5.0e4 = 0.5 + 0.5 + 0.6 = 1.6 m^3
  // The decoy cells would contribute ~1e10 if they were included.
  CHECK_CLOSE(1.6, elmColumnWaterVolume(surf_wc, surf_n_liq, wc, n_liq, cells), TIGHT);
}

// The surface cell must use the *surface* molar density, not the subsurface
// one.  Identical mol amounts with different densities must not cancel.
void
ColumnVolumeSurfaceUsesSurfaceDensity()
{
  const double wc[1] = { 1.0e4 };      // mol
  const double n_liq[1] = { 1.0e4 };   // mol m-3 -> 1.0 m^3
  const std::vector<int> cells = { 0 };
  const double surf_wc = 1.0e4;        // same mol amount...
  const double surf_n_liq = 2.0e4;     // ...but half the volume: 0.5 m^3

  CHECK_CLOSE(1.5, elmColumnWaterVolume(surf_wc, surf_n_liq, wc, n_liq, cells), TIGHT);
}

// Storage depth divides by surface area, so a non-unit area is exercised:
// 2.5 m^3 * 1000 kg/m^3 / 250 m^2 = 10 mm.
void
StorageDepthNonUnitArea()
{
  CHECK_CLOSE(10.0, elmStorageDepth(2.5, 250.0), TIGHT);
}

// Full driver chain on a synthetic balanced column: storage rises by exactly
// the source input over the step, so errh2o must be ~0.  This ties the mol->mm
// conversion to the m/s->mm/s flux scaling: dropping the 1e3 factor, or the
// denh2o/area division, breaks the cancellation.
void
DriverChainBalancedColumn()
{
  const double dt = 1800.0;            // s
  const double area = 250.0;           // m^2
  const double n_liq[2] = { 5.0e4, 5.0e4 };
  const double surf_n_liq = 5.0e4;
  const std::vector<int> cells = { 0, 1 };

  // Begin: 2.0 m^3 subsurface (1.0 each), no surface water.
  const double wc_old[2] = { 5.0e4, 5.0e4 };
  const double swc_old = 0.0;
  const double begwb = elmStorageDepth(
    elmColumnWaterVolume(swc_old, surf_n_liq, wc_old, n_liq, cells), area);
  CHECK_CLOSE(8.0, begwb, TIGHT);      // 2.0 * 1000 / 250

  // A source of 1e-6 m/s over 1800 s adds 1.8e-3 m depth = 1.8 mm, which over
  // 250 m^2 is 0.45 m^3 = 2.25e4 mol.  Put it in the surface cell.
  const double source_mps = 1.0e-6;
  const double swc_new = 2.25e4;
  const double endwb = elmStorageDepth(
    elmColumnWaterVolume(swc_new, surf_n_liq, wc_old, n_liq, cells), area);
  CHECK_CLOSE(9.8, endwb, TIGHT);      // 8.0 + 1.8

  const double errh2o = elmWaterBalanceError(
    endwb, begwb, source_mps * ELM_M_PER_S_TO_MM_PER_S, 0.0, 0.0, 0.0, 0.0, dt);
  CHECK_CLOSE(0.0, errh2o, 1.0e-10);
}

// Same chain, but water leaves as evaporation while storage is unchanged, so
// the imbalance is exactly the sink over the step.
//
// NOTE: this pins the sign convention as currently coded -- sinks
// (evap/tran/baseflow/runoff) are positive-leaving, matching the *_mps coupling
// convention.  That convention is still flagged for runtime validation against
// a real run (see the NOTE on checkELMWaterBalance_); this case exists so an
// accidental sign flip in the conversion fails a test, not to assert that the
// assumption itself has been confirmed.
void
DriverChainEvaporativeLeak()
{
  const double dt = 1800.0;
  const double area = 250.0;
  const double n_liq[1] = { 5.0e4 };
  const double surf_n_liq = 5.0e4;
  const std::vector<int> cells = { 0 };
  const double wc[1] = { 5.0e4 };      // 1.0 m^3, unchanged over the step

  const double wb = elmStorageDepth(
    elmColumnWaterVolume(0.0, surf_n_liq, wc, n_liq, cells), area);
  CHECK_CLOSE(4.0, wb, TIGHT);         // 1.0 * 1000 / 250

  // Storage did not change, so all of the evaporated water is unaccounted for:
  // errh2o = 0 - (-evap)*dt = +evap*dt = 1e-6 * 1e3 * 1800 = 1.8 mm.
  const double evap_mps = 1.0e-6;
  const double errh2o = elmWaterBalanceError(
    wb, wb, 0.0, evap_mps * ELM_M_PER_S_TO_MM_PER_S, 0.0, 0.0, 0.0, dt);
  CHECK_CLOSE(1.8, errh2o, 1.0e-10);
}

struct TestCase {
  const char* name;
  void (*fn)();
};

const TestCase TESTS[] = {
  { "PerfectBalance", PerfectBalance },
  { "PureSource", PureSource },
  { "EvapSink", EvapSink },
  { "TranSink", TranSink },
  { "RunoffSink", RunoffSink },
  { "BaseflowSink", BaseflowSink },
  { "EvapPlusTran", EvapPlusTran },
  { "StorageOnlyZeroDt", StorageOnlyZeroDt },
  { "CombinedNonZero", CombinedNonZero },
  { "ColumnVolumeNonUniformDensity", ColumnVolumeNonUniformDensity },
  { "ColumnVolumeScatteredCellIds", ColumnVolumeScatteredCellIds },
  { "ColumnVolumeSurfaceUsesSurfaceDensity", ColumnVolumeSurfaceUsesSurfaceDensity },
  { "StorageDepthNonUnitArea", StorageDepthNonUnitArea },
  { "DriverChainBalancedColumn", DriverChainBalancedColumn },
  { "DriverChainEvaporativeLeak", DriverChainEvaporativeLeak },
};

} // namespace

int
main(int /*argc*/, char* /*argv*/[])
{
  int failed_tests = 0;
  for (const TestCase& t : TESTS) {
    const int before = g_failures;
    t.fn();
    if (g_failures > before) {
      ++failed_tests;
      std::fprintf(stderr, "[FAIL] %s\n", t.name);
    } else {
      std::fprintf(stdout, "[PASS] %s\n", t.name);
    }
  }

  const int ntests = static_cast<int>(sizeof(TESTS) / sizeof(TESTS[0]));
  std::fprintf(stdout,
               "ELM_ATS_WATER_BALANCE_KERNEL: %d/%d tests passed "
               "(%d checks, %d failures)\n",
               ntests - failed_tests, ntests, g_checks, g_failures);

  return g_failures == 0 ? 0 : 1;
}
