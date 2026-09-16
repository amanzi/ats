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
