/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

A constitutive relation for the viscosity of water as a function of temperature
in K, given as an empirical series expansion fit to data.

Used by setting

`"viscosity relation type`" = `"liquid water`"

.. _viscosity-water-spec:
.. admonition:: viscosity-water-spec

   NONE

*/

#ifndef AMANZI_RELATIONS_VISCOSITY_WATER_HH_
#define AMANZI_RELATIONS_VISCOSITY_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "dbc.hh"
#include "viscosity_relation.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class ViscosityWater : public ViscosityRelation {
 public:
  explicit ViscosityWater(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T);
  virtual double DViscosityDT(double T);

 protected:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of viscosity < T1
  const double kav1_, kbv1_, kcv1_;

  // -- temperature dependence of viscosity > T1
  const double kbv2_, kcv2_, kT1_;

 private:
  static Utils::RegisteredFactory<ViscosityRelation, ViscosityWater> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
