/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

*/

#ifndef AMANZI_RELATIONS_EOS_WATER_HH_
#define AMANZI_RELATIONS_EOS_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSWater : public EOSConstantMolarMass {
 public:
  explicit EOSWater(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(std::vector<double>& params) override;
  virtual double DMassDensityDT(std::vector<double>& params) override;
  virtual double DMassDensityDp(std::vector<double>& params) override;
  virtual double DMassDensityDC(std::vector<double>& params) override { return 0; }

  virtual bool IsConcentration() override { return false; }
  virtual bool IsTemperature() override { return true; }
  virtual bool IsPressure() override { return true; }

 private:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_, kd_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactory<EOS, EOSWater> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
