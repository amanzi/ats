/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

*/

#ifndef AMANZI_RELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_
#define AMANZI_RELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_


#include "Teuchos_ParameterList.hpp"
#include "Factory.hh"
#include "vapor_pressure_relation.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class VaporPressureWater : public VaporPressureRelation {
 public:
  explicit VaporPressureWater(Teuchos::ParameterList& plist);

  virtual double SaturatedVaporPressure(double T);
  virtual double DSaturatedVaporPressureDT(double T);

 private:
  Teuchos::ParameterList plist_;
  const double ka0_;
  const double ka_, kb_, kc_, kd_;

  static Utils::RegisteredFactory<VaporPressureRelation, VaporPressureWater> factory_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
