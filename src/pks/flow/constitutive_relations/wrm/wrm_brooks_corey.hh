/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef FLOWRELATIONS_WRM_BROOKS_COREY_
#define FLOWRELATIONS_WRM_BROOKS_COREY_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMBrooksCorey : public WRM {
 public:
  explicit WRMBrooksCorey(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double pc);
  double d_k_relative(double pc);

  double saturation(double pc);
  double d_saturation(double pc);

  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

 private:
  void InitializeFromPlist_();

 private:

  Teuchos::ParameterList& plist_;

  double lambda_, l_, alpha_;  // Brooks and Corey parameters: lambda, alpha
  double sr_;  // residual saturation
  int function_;

  double pc0_;  // regularization threshold (usually 0 to 500 Pa)
  double a_, b_, factor_, pc_bubble_;  // frequently used constant

  static Utils::RegisteredFactory<WRM,WRMBrooksCorey> factory_;

};

}  // namespace
}  // namespace
}  // namespace


#endif
