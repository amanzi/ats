/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Basic land cover/plant function type
/*!

This is a simple base class for PFTs, but currently called LandCover_t to
differentiate from PFT, which is currently used by the dynamic BGC PK.

This class gets used for compositional land cover, such as NLCD data.  To use
this, a list of these specs is included in the "initial conditions" list of
State, and can be used by multiple evaluators somehow.  These are much cleaner
to work with in "new state."

The names of each spec in the list correspond to the region on the surface mesh
of that land cover type.  These regions should partition the surface domain,
although this is not checked, and exact behavior depends upon the evaluators
using these.  In practice, the regions are processed in the order they appear
in the list, and if a cell is in multiple regions, the last region's value will
take precedence.

In general, all parameter values are initialized as signalling NaNs, so if they
are used but not set, the code with throw.  This allows users to _not_ include
values if they aren't using models that need those values.

The parameters are used in a few different models.  In general, the use of this
class allows multiple evaluators to share these parameters.  The advantage of
this approach is that there is no fear of inconsistent parameters across
models.  The downside of this approach is that these parameters must have the
same region-based partitioning.

.. _land-cover-spec:
.. admonition:: land-cover-spec

   * `"rooting depth max [m]`" ``[double]`` **NaN** Below this the rooting
     fraction is set to 0. [m]
   * `"rooting profile alpha [-]`" ``[double]`` **NaN** alpha in the rooting
     profile function [-]
   * `"rooting profile beta [-]`" ``[double]`` **NaN** beta in the rooting
     profile function [-] Note that these are from the CLM 4.5 Technical Note.

   * `"capillary pressure at fully closed stomata [Pa]`" ``[double]`` **NaN**
   * `"capillary pressure at fully open stomata [Pa]`" ``[double]`` **NaN**
     Transpiration is typically downregulated by a limiter that is empirically
     modeling stomata closure.  Typically it varies linearly from 0 to 1 as a
     function of capillary pressure, between these two values.  Note that
     these should be positive! [Pa]

   * `"maximum xylem capillary pressure [Pa]`" ``[double]`` **NaN**
     Maximum capillary pressure at the plant collar before the plant starts to close stomata.

   * `"leaf on time [doy]`" ``[double]`` **NaN** Day of year, relative to time
     0, when leaves begin transpiring.  Note that -1 implies evergreen. [doy]
   * `"leaf off time [doy]`" ``[double]`` **NaN** Day of year, relative to
     time 0, when leaves stop transpiring.  Note that -1 implies
     evergreen. [doy]

     Note that `"leaf on doy`" and `"leaf off doy`" are relative to the
     simulation's zero time, not the start time.  Typically these are Julian
     day of the year, but this assumes that the 0 time of the simulation (not
     the "start time", but time 0!) is Jan 1.  This leaf on/off cycle is modulo
     the `"year duration`" (typically 1 noleap).  Note if `"leaf off doy`" <
     `"leaf on time`" is ok too -- this is the case if simulation time zero is
     mid-summer.

   * `"Priestley-Taylor alpha of snow [-]`" ``[double]`` **NaN** Evaporation
     coefficient in the Priestley-Taylor model, used in sublimation of snow.
   * `"Priestley-Taylor alpha of bare ground [-]`" ``[double]`` **NaN**
     Evaporation coefficient in the Priestley-Taylor model, used in bare soil
   * `"Priestley-Taylor alpha of canopy [-]`" ``[double]`` **NaN** Evaporation
     coefficient in the Priestley-Taylor model, used in sublimation of snow.
   * `"Priestley-Taylor alpha of transpiration [-]`" ``[double]`` **NaN**
     Evaporation coefficient in the Priestley-Taylor model, used in
     transpiration.

   * `"interception coefficient [-]`" ``[double]`` **NaN** Fraction, per unit
     LAI, of water intercepted by the canopy.

   * `"albedo of bare ground [-]`" ``[double]`` **NaN** Albedo of the land
     cover type, ranging from [0,1]
   * `"emissivity of bare ground [-]`" ``[double]`` **NaN** Emissivity of the
     land cover type's bare ground, ranging from [0,1]
   * `"albedo of canopy [-]`" ``[double]`` **NaN** Albedo of the land cover
     type's canopy, ranging from [0,1]

   * `"Beer's law extinction coefficient, shortwave [-]`"  ``[double]``  **NaN**
   * `"Beer's law extinction coefficient, longwave [-]`" ``[double]`` **NaN**
     Beer's law provides the attenuation of light through a diffuse medium (in
     this case leaves) as a function of concentration/density (in this case
     LAI) and an extinction coefficient, sometimes called the absorptivity,
     which is usually a function of wavelength, but here is split into long
     and shortwave constants.

   * `"snow transition depth [m]`" ``[double]`` **NaN** See below.
   * `"water transition depth [m]`" ``[double]`` **NaN** Due to the need for
     smooth transitions (numerically) between various surface balance
     conditions, frequently a subgrid-scale model is used that, rather than
     trying to smooth each model, varies an area fraction between 1 and 0 for
     multiple components (e.g. snow-covered, water-covered, bare-ground) of
     the surface cell.  These area fractions are determined by a model, but
     most frequently are a linear variation from 1 (e.g. fully covered by
     water) at this ponded depth and above, and 0 (e.g. no area covered by
     water) at 0 ponded depth.  Then, multiple surface energy balances (SEBs)
     are calculated, one assuming the ponded depth is this thick, and the
     various components of the SEB are averaged by that area fraction.

   * `"roughness length of bare ground [m]`" ``[double]`` **NaN** Roughness
     length of the bare soil, used in calculating sensible/latent heat in the
     physics-based SEB model.  A typical value is 0.04.
   * `"roughness length of snow [m]`" ``[double]`` **NaN** Roughness length of
     the snow-covered soil, used in calculating sensible/latent heat in the
     physics-based SEB model.  A typical value is 0.004.


   - `"Manning's n [?]`" ``[double]`` **NaN** Manning's n [??]  THIS IS NOT
     CURRENTLY USED.

*/

#pragma once

#include <map>
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace SurfaceBalance {

//
// LandCover, or PFT struct
//
struct LandCover {
  LandCover(Teuchos::ParameterList& plist);

  // rooting profiles
  double rooting_depth_max;
  double rooting_profile_alpha;
  double rooting_profile_beta;

  // parameters in the transpiration reduction function
  double stomata_closed_capillary_pressure;
  double stomata_open_capillary_pressure;

  double maximum_xylem_capillary_pressure;

  // manning's coef
  double mannings_n;

  // transpiration controls
  double leaf_on_doy;
  double leaf_off_doy;

  // priestley-taylor model parameters
  double pt_alpha_snow;
  double pt_alpha_canopy;
  double pt_alpha_ground;
  double pt_alpha_transpiration;

  // radiation parameters
  double albedo_ground;
  double albedo_canopy;
  double emissivity_ground;

  double beers_k_sw;
  double beers_k_lw;

  // transition thickness between snow and bare ground
  // likely this is a property of the understory veg
  double snow_transition_depth;  // [m]
  double water_transition_depth; // [m]

  // soil properties controlling evaporation
  double roughness_ground; // [m] Fetch length for latent/sensible heat fluxes.
  double roughness_snow;   // [m] Fetch length for latent/sensible heat fluxes.
};

// this one includes error checking for NaNs
using LandCoverMap = std::map<std::string, LandCover>;
LandCoverMap getLandCover(Teuchos::ParameterList plist,
                          const std::vector<std::string>& required_pars);

namespace Impl {

void checkValid(const std::string& region, const LandCover& lc, const std::string& parname);
LandCoverMap getLandCover(Teuchos::ParameterList& plist);


} // namespace Impl


} // namespace SurfaceBalance
} // namespace Amanzi
