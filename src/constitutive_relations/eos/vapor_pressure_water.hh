/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_
#define AMANZI_RELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_


#include "Teuchos_ParameterList.hpp"
#include "Factory.hh"
#include "vapor_pressure_relation.hh"

namespace Amanzi {
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
} // namespace Amanzi

#endif
