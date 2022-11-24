/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Basic interface of a Viscosity.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_VISCOSITY_MODEL_HH_
#define AMANZI_RELATIONS_VISCOSITY_MODEL_HH_

namespace Amanzi {
namespace Relations {

// Equation of State model
class ViscosityRelation {
 public:
  virtual ~ViscosityRelation(){};

  // Virtual methods that form the Viscosity
  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;
};

} // namespace Relations
} // namespace Amanzi

#endif
