/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for TC implementations.
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_factory.hh"

// explicity instantitate the static data of Factory<EOS>
template<>
Amanzi::Utils::Factory<Amanzi::ATS_Physics::Energy::ThermalConductivityThreePhase>::map_type*
  Amanzi::Utils::Factory<Amanzi::ATS_Physics::Energy::ThermalConductivityThreePhase>::map_;
