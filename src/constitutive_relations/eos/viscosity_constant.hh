/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Constant prescribed viscosity.

  http://software.lanl.gov/ats/trac

*/


#ifndef AMANZI_RELATIONS_VISCOSITY_CONSTANT_HH_
#define AMANZI_RELATIONS_VISCOSITY_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Factory.hh"
#include "dbc.hh"
#include "viscosity_relation.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// Equation of State model
class ViscosityConstant : public ViscosityRelation {
 public:
  explicit ViscosityConstant(Teuchos::ParameterList& visc_plist);

  virtual double Viscosity(double T) { return visc_; }
  virtual double DViscosityDT(double T) { return 0.; }

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList visc_plist_;
  double visc_;

 private:
  static Utils::RegisteredFactory<ViscosityRelation, ViscosityConstant> factory_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
