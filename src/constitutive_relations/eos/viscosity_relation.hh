/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Basic interface of a Viscosity Relation.

*/

#ifndef AMANZI_RELATIONS_VISCOSITY_MODEL_HH_
#define AMANZI_RELATIONS_VISCOSITY_MODEL_HH_

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// Equation of State model
class ViscosityRelation {
 public:
  virtual ~ViscosityRelation() {};

  // Virtual methods that form the Viscosity
  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
