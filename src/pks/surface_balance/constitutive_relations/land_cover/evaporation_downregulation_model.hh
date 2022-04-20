/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Downregulates bare soil evaporation through a dessicated zone.
/*!

Calculates evaporative resistance through a dessicated zone.

Sakagucki and Zeng 2009 equations 9 and 10.

Requires two parameters,

* `"dessicated zone thickness [m]`" Thickness over which vapor must diffuse
  when the soil is dry.

* `"Clapp and Hornberger b [-]`" Exponent of the Clapp & Hornberger curve for
  the top layer of soil.  Nominally this could probably be pulled from van
  Genuchten curves that we typically use, but it doesn't appear to be the most
  important parameter.

These may be provided via parameter list or LandCover type.

*/

#pragma once

#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporationDownregulationModel {

 public:
  explicit
  EvaporationDownregulationModel(double dessicated_zone_thickness,
          double clapp_horn_b);
  EvaporationDownregulationModel(Teuchos::ParameterList& plist);
  EvaporationDownregulationModel(const LandCover& lc);

  double Evaporation(double sg, double poro, double pot_evap) const;

  double DEvaporationDSaturationGas(double sg, double poro, double pot_evap) const;
  double DEvaporationDPorosity(double sg, double poro, double pot_evap) const;
  double DEvaporationDPotentialEvaporation(double sg, double poro, double pot_evap) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double dess_dz_;
  double Clapp_Horn_b_;

};

} //namespace
} //namespace
} //namespace

