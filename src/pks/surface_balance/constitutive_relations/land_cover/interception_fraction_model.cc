/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Fraction of incoming water that is intercepted.

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "errors.hh"
#include "interception_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
InterceptionFractionModel::InterceptionFractionModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
InterceptionFractionModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  alpha_ = plist.get<double>("leaf area interception fraction [-]", 0.25);
  if (alpha_ < 0 && alpha_ > 1) {
    Errors::Message msg;
    msg << "InterceptionFraction: invalid \"leaf area interception fraction [-]\", must be in [0,1] (provided: "
        << alpha_ << ")";
    Exceptions::amanzi_throw(alpha_);
  }
}


// main method
double
InterceptionFractionModel::InterceptionFraction(double ai) const
{
  return alpha_*(1 - exp(-0.5*ai));
}

double
InterceptionFractionModel::DInterceptionFractionDAreaIndex(double ai) const
{
  return 0.5*alpha_*exp(-0.5*ai);
}

} //namespace
} //namespace
} //namespace

