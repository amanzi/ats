/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  WRM which calls another WRM for saturation but sets 0 rel perm.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "dbc.hh"
#include "wrm_linear_relperm.hh"

namespace Amanzi {
namespace Flow {

WRMLinearRelPerm::WRMLinearRelPerm(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
};


void
WRMLinearRelPerm::InitializeFromPlist_()
{
  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");

  WRMFactory fac;
  wrm_ = fac.createWRM(sublist);
};


} // namespace Flow
} // namespace Amanzi
