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

#ifndef _FLOWRELATIONS_WRM_LINEAR_SYSTEM_
#define _FLOWRELATIONS_WRM_LINEAR_SYSTEM_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMLinearSystem : public WRM {

public:
  explicit WRMLinearSystem(Teuchos::ParameterList& plist);

  // required methods from the base class
  virtual double k_relative(double pc) { return 1.0; }
  virtual double d_k_relative(double pc) { return 0.0; }
  virtual double saturation(double pc) { return sat_at_zero_pc_ + alpha_*pc; }
  virtual double d_saturation(double pc) { return alpha_; }
  virtual double capillaryPressure(double saturation) { return (saturation-sat_at_zero_pc_)/alpha_; }
  virtual double d_capillaryPressure(double saturation) { return 1./alpha_; }
  virtual double residualSaturation() { return 0.0; }

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;
  double alpha_;
  double sat_at_zero_pc_;

  static Utils::RegisteredFactory<WRM,WRMLinearSystem> factory_;
};

} //namespace
} //namespace

#endif
