/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rich Fiorella, Ethan Coon
*/

//! Unit test for TimeAdvancer step rejection/retry driven by the ELM water balance.
/*!

  Exercises the reject-and-retry branch of TimeAdvancer::advance() that is
  taken when the optional step validity check (set by ELM_ATSDriver to gate
  step acceptance on ELM's water mass balance) returns false for a step the
  PK otherwise converged.

  No mesh, physics PK, or input file is needed.  A stub PK always converges
  with a fixed requested dt and records every Advance/Commit/Fail call.  The
  validity check calls the real ELM kernel elm_ats_water_balance_error_c with
  a synthetic constant leak L (no storage change), so errh2o = L * dt grows
  with dt.  With a tolerance between L*dt and L*dt/2 the first step must be
  rejected, dt reduced, and the retry accepted.

*/

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>
#include <utility>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

#include "VerboseObject_objs.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "TimeStepManager.hh"

#include "time_advancer.hh"
#include "elm_balance_interface_private.hh"

namespace {

// --- Minimal dependency-free test harness ----------------------------------

int g_checks = 0;
int g_failures = 0;

void
check(bool ok, const char* file, int line, const char* expr)
{
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::fprintf(stderr, "%s:%d: FAILED CHECK(%s)\n", file, line, expr);
  }
}

#define CHECK(expr) check((expr), __FILE__, __LINE__, #expr)

bool
close(double a, double b)
{
  return std::fabs(a - b) <= 1.0e-10 * (1.0 + std::fabs(a) + std::fabs(b));
}

using Interval = std::pair<double, double>;

bool
same(const std::vector<Interval>& actual, const std::vector<Interval>& expected)
{
  if (actual.size() != expected.size()) return false;
  for (std::size_t i = 0; i != actual.size(); ++i) {
    if (!close(actual[i].first, expected[i].first) || !close(actual[i].second, expected[i].second))
      return false;
  }
  return true;
}

// --- Stub PK ----------------------------------------------------------------

// Always converges with a fixed requested dt; records call history.
class StubPK : public Amanzi::PK {
 public:
  StubPK(Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
         const Teuchos::RCP<Amanzi::State>& S,
         double dt)
    : Amanzi::PK(pk_tree, global_plist, S, Teuchos::null), dt_(dt)
  {}

  void parseParameterList() override {}
  void Setup() override {}
  void Initialize() override {}

  double get_dt() override { return dt_; }
  void set_dt(double dt) override { dt_ = dt; }

  bool AdvanceStep(double t_old, double t_new, bool reinit) override
  {
    advanced.emplace_back(t_old, t_new);
    return false;
  }
  void CommitStep(double t_old, double t_new, const Amanzi::Tag& tag) override
  {
    committed.emplace_back(t_old, t_new);
  }
  void FailStep(double t_old, double t_new, const Amanzi::Tag& tag) override
  {
    failed.emplace_back(t_old, t_new);
  }

  void State_to_Solution(const Amanzi::Tag& tag, Amanzi::TreeVector& soln) override {}
  void Solution_to_State(const Amanzi::TreeVector& soln, const Amanzi::Tag& tag) override {}

  std::vector<Interval> advanced, committed, failed;

 private:
  double dt_;
};

// --- Fixture ----------------------------------------------------------------

constexpr double PK_DT = 100.0;   // [s] dt the PK always requests
constexpr double LEAK = 1.0e-3;   // [mm s-1] synthetic imbalance
constexpr double T_END = 100.0;   // [s]

struct RunResult {
  std::vector<double> checked_dts;  // dt seen by each validity-check call
  std::vector<Interval> advanced, committed, failed;
  double t_final;
  int cycles;
};

// Runs one outer advance over [0, T_END] with the ELM-kernel validity check at
// tolerance tol [mm].  reduction <= 0 keeps the TimeAdvancer default.
RunResult
run(double tol, double reduction = -1.0)
{
  using Amanzi::Tags::CURRENT;
  using Amanzi::Tags::NEXT;

  auto global_plist = Teuchos::rcp(new Teuchos::ParameterList("main"));
  global_plist->sublist("PKs").sublist("stub");
  Teuchos::ParameterList pk_tree("stub");

  Teuchos::ParameterList state_plist("state");
  auto S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->require_time(CURRENT);
  S->require_time(NEXT);

  auto pk = Teuchos::rcp(new StubPK(pk_tree, global_plist, S, PK_DT));
  auto tsm = Teuchos::rcp(new Amanzi::Utils::TimeStepManager());
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("TimeAdvancerRetryTest", "low"));

  auto ta_plist = Teuchos::rcp(new Teuchos::ParameterList("cycle driver"));
  if (reduction > 0.0) ta_plist->set("validity timestep reduction factor", reduction);

  ATS::TimeAdvancer ta(ta_plist, S, pk, tsm, CURRENT, NEXT, vo, Teuchos::null);

  RunResult r;
  ta.set_step_validity_check([&r, tol](double t_old, double t_new) {
    const double dt = t_new - t_old;
    r.checked_dts.push_back(dt);
    const double begwb = 100.0, endwb = begwb, zero = 0.0, leak = LEAK;
    const double errh2o =
      elm_ats_water_balance_error_c(&endwb, &begwb, &zero, &leak, &zero, &zero, &zero, &dt);
    return std::abs(errh2o) <= tol;
  });

  ta.setup();
  S->Setup();
  S->set_time(Amanzi::Tags::DEFAULT, 0.0);
  S->set_time(CURRENT, 0.0);
  S->set_time(NEXT, 0.0);
  S->set_cycle(Amanzi::Tags::DEFAULT, 0);
  S->set_cycle(NEXT, 0);
  ta.initialize();

  ta.advance(0.0, T_END);

  r.advanced = pk->advanced;
  r.committed = pk->committed;
  r.failed = pk->failed;
  r.t_final = S->get_time(CURRENT);
  r.cycles = S->get_cycle(NEXT);
  return r;
}

// --- Test cases -------------------------------------------------------------

// errh2o(dt=100) = 0.1 mm > tol, errh2o(dt=50) = 0.05 mm <= tol: the first
// step is rejected, retried at half dt, then the remainder is taken.
void
RejectThenRetry()
{
  const RunResult r = run(0.075);
  CHECK(r.checked_dts.size() == 3);
  if (r.checked_dts.size() == 3) {
    CHECK(close(r.checked_dts[0], 100.0));
    CHECK(close(r.checked_dts[1], 50.0));
    CHECK(close(r.checked_dts[2], 50.0));
  }
  CHECK(same(r.advanced, { { 0., 100. }, { 0., 50. }, { 50., 100. } }));
  CHECK(same(r.failed, { { 0., 100. } }));
  CHECK(same(r.committed, { { 0., 50. }, { 50., 100. } }));
  CHECK(close(r.t_final, T_END));
  CHECK(r.cycles == 2);
}

// The "validity timestep reduction factor" parameter sets the retry dt.
void
ReductionFactor()
{
  const RunResult r = run(0.075, 0.25);
  CHECK(r.checked_dts.size() >= 2);
  if (r.checked_dts.size() >= 2) {
    CHECK(close(r.checked_dts[0], 100.0));
    CHECK(close(r.checked_dts[1], 25.0));
  }
  CHECK(r.failed.size() == 1);
  CHECK(close(r.t_final, T_END));
}

// Control: a loose tolerance accepts the first step, so any rejection above
// comes from the validity check.
void
LooseToleranceAccepts()
{
  const RunResult r = run(1.0);
  CHECK(r.checked_dts.size() == 1);
  CHECK(r.failed.empty());
  CHECK(same(r.committed, { { 0., 100. } }));
  CHECK(close(r.t_final, T_END));
  CHECK(r.cycles == 1);
}

struct TestCase {
  const char* name;
  void (*fn)();
};

const TestCase TESTS[] = {
  { "RejectThenRetry", RejectThenRetry },
  { "ReductionFactor", ReductionFactor },
  { "LooseToleranceAccepts", LooseToleranceAccepts },
};

} // namespace

int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  int failed_tests = 0;
  {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    for (const TestCase& t : TESTS) {
      const int before = g_failures;
      try {
        t.fn();
      } catch (const std::exception& e) {
        ++g_failures;
        std::fprintf(stderr, "%s threw: %s\n", t.name, e.what());
      }
      if (g_failures > before) {
        ++failed_tests;
        std::fprintf(stderr, "[FAIL] %s\n", t.name);
      } else {
        std::fprintf(stdout, "[PASS] %s\n", t.name);
      }
    }

    const int ntests = static_cast<int>(sizeof(TESTS) / sizeof(TESTS[0]));
    std::fprintf(stdout,
                 "ELM_ATS_TIME_ADVANCER_VALIDITY_RETRY: %d/%d tests passed "
                 "(%d checks, %d failures)\n",
                 ntests - failed_tests, ntests, g_checks, g_failures);
  }
  Kokkos::finalize();
  return g_failures == 0 ? 0 : 1;
}
