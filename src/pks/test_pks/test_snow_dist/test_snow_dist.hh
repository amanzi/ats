/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Test PK for testing snow distribution

------------------------------------------------------------------------- */

#ifndef PKS_TEST_SNOW_DIST_HH_
#define PKS_TEST_SNOW_DIST_HH_

#include "pk_factory_ats.hh"
#include "pk_physical_base.hh"

namespace Amanzi {

class TestSnowDist : public PKPhysicalBase {
 public:
  TestSnowDist(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               Teuchos::ParameterList& FElist,
               const Teuchos::RCP<TreeVector>& solution)
    : PKDefaultBase(plist, FElist, solution), PKPhysicalBase(plist, FElist, solution)
  {
    plist_->set("solution key", "snow_depth");
    sink_type_ = plist_->get<std::string>("sink type", "none");
    sink_value_ = plist_->get<double>("sink value", 1.);
  }

  // ConstantTemperature is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool advance(double dt);

  double get_dt()
  {
    return 24 * 60 * 60.0; // 1 day
  }

 protected:
  std::string sink_type_;
  double sink_value_;

 private:
  // factory registration
  static RegisteredPKFactory_ATS<TestSnowDist> reg_;
};

} // namespace Amanzi

#endif
