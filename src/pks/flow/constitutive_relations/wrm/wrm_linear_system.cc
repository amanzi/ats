/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

*/
//! A linear sat-pc curve.

/*!

  A linear sat-pc curve, plus a constant rel perm, makes the system linear, so
  nonlinear solver should always converge in one step.

  No error-checking, so the user is responsible for ensuring that the pressure
  is always less than atmospheric and within the acceptable range of the slope.

  Note this is mostly for testing.

*/

#include "wrm_linear_system.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMLinearSystem::WRMLinearSystem(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

void WRMLinearSystem::InitializeFromPlist_() {
  sat_at_zero_pc_ = plist_.get<double>("saturation at pc=0", 1.0);
  if (plist_.isParameter("alpha")) {
    alpha_ = plist_.get<double>("alpha");
  } else {
    double max_pc = plist_.get<double>("max pc");
    alpha_ = -sat_at_zero_pc_/max_pc;
  }
};

}  // namespace
}  // namespace
