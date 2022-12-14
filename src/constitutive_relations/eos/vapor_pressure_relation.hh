/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas air with a molar fraction of water vapor.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_VAPOR_PRESSURE_RELATION_HH_
#define AMANZI_RELATIONS_EOS_VAPOR_PRESSURE_RELATION_HH_

namespace Amanzi {
namespace Relations {

class VaporPressureRelation {
 public:
  virtual ~VaporPressureRelation(){};

  virtual double SaturatedVaporPressure(double T) = 0;
  virtual double DSaturatedVaporPressureDT(double T) = 0;
};

} // namespace Relations
} // namespace Amanzi

#endif
