/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*----------------------------------------------------------------------------

ATS

Functions needed for using CLM as a land surface model.

Usage:
=========

Initialization:
    (order matters!)
-----------------

init()
set_zero_time()
set_initial_state()
set_ground_properties()
setup_begin()
set_dz()
set_et_controls()
set_irrigation_controls()       [optional, not currently implemented FIXME]
setup_end()
set_dz()                        [a second time, FIXME]

Advance timestep:
   (advance_time() last!)
-------------------
set_tksat_from_porosity()       [if porosity changes, alternatively call
                                 once after initialization]
set_wc()
set_pressure()                  [can these be merged? FIXME]
set_met_data()                  [figure out expected Met data! FIXME]
advance_time()


Forcing for hydrologic model:
-------------------------------
get_total_mass_fluxes()


Diagnostics:
--------------
get_total_energy_fluxes()
get_mass_fluxes()
get_diagnostics()



----------------------------------------------------------------------------*/


#ifndef ATS_CLM_INTERFACE_HH_
#define ATS_CLM_INTERFACE_HH_

#include <cstdint>
#include <vector>

#include "Epetra_MultiVector.h"

#define NUM_LC_CLASSES 18

namespace ATS {
namespace CLM {

//
// Begin initialization, allocating space for driver, grid.
// ------------------------------------------------------------------
// Input:
//   ncells     | number of grid cells in the subsurface
//   ncolumns   | number of columns.  Note it is currently assumed that
//              |  ncolumns divides ncells evenly, with a fixed number of
//              |  cells per column.  (FIXME)
//   startcode  | 1 if new run, 2 if restart (not currently supported?)
//   rank       | rank (for logging purposes only)
//   verbosity  | 0 (None) - 1 (Low) - 2 (High) - 3 (Extreme)
//
int init(int ncells, int ncolumns, int startcode, int rank, int verbosity);


//
// Set the reference, time 0, in years.
// ------------------------------------------------------------------
// Input:
//   zero_time  | The time from which all time is measured. [year]
//              |  i.e. 2010.0 for midnight Jan 1, 2010
//              |  NOTE: what does CLM actually need this for?
//              |  CLM calculates zeniths and incident radiation somehow,
//              |  but I would have thought this was assumed to be given
//              |  by incoming radiation.  So what incoming met data
//              |  temporal resolution does CLM expect? --etc
//
int set_zero_time(double zero_time);

//
// Sets the initial, assumed uniform, state.
// ------------------------------------------------------------------
// Input:
//   temperature| Soil, snow, and water temperature [K]
//   snow depth | Initial snow depth [m]
//
int set_initial_state(double temperature, double snow_depth);


//
// Set soil properties
// ------------------------------------------------------------------
// Input:
//   latlon     | lat/lon.  +Northern and Eastern hemispheres. [degrees]
//              |  Array of shape [ncolumns, 2]
//   sand,clay  | soil texture, fractions (must range from 0-1) [-]
//              |  Size ncells.
//   color_index| index into soil color models?  See clm1d_varcon albsat
//              |  and albdry.  Size ncolumns.
//   f_ground   | Fraction of each land type.  Land types are set in drv_vegp.dat
//              |  Array of shape [ncolumns, NUM_LC_CLASSES]
//              |  Note only the dominant land type is used, as maxt is
//              |  hard-coded 1 and this assumption is used in a few places in drv.
//
int set_ground_properties(double* latlon,
                          const Epetra_MultiVector& sand,
                          const Epetra_MultiVector& clay,
                          const std::vector<int>& color_index,
                          double* fractional_ground);

//
// Begin setup. Sets the tiles and pushes info from driver into grid/tile
// ------------------------------------------------------------------
//
int setup_begin();


//
// sets the cell thicknesses in the vertical
// ------------------------------------------------------------------
// Input:
//   dz         | vector of cell dz [m]  Size ncells.
//
int set_dz(const Epetra_MultiVector& dz);


//
// Sets ET controls
// ------------------------------------------------------------------
// Input:
//   beta_type          | ??
//   veg_water_stress_type
//                      | ??
//   wilting_point      | ??
//   field_capacity     | ??
//   res_sat            | residual saturation [-]
//                      |  make this variable across columns?  cells? FIXME
//
int set_et_controls(int beta_type,
                    int veg_water_stress_type,
                    double wilting_point,
                    double field_capacity,
                    double res_sat);

//
// End setup. Pushes grid, tile, drv info into clm1d column instances.
// ------------------------------------------------------------------
//
int setup_end();


//
// Move water content data into clm1d columns
// ------------------------------------------------------------------
// Input:
//   porosity   | [-]  Size ncells.
//   saturation | [-]  Size ncells.
//
int set_wc(const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation);


//
// Set saturated thermal conductivity from porosity.
// -----------------------------------------------------------------------
//  NOTE: This is seperate from set_wc() to allow it to be called once and
//   assume to be fixed throughout the simulation (small porosity
//   variation due to compressibility only.
//
// Input:
//   porosity   | [-]  Size ncells.
//
int set_tksat_from_porosity(const Epetra_MultiVector& porosity);


//
// Move pressure data into 1d columns.
// ------------------------------------------------------------------
// Input:
//   pressure   | [Pa] Size ncells.
//   p_atm      | Atmospheric pressure.  [Pa] Note CLM works in head
//              |  units, so p_atm is used to convert Pa to mm.  It
//              |  isn't obvious that this is a good idea in cases of
//              |  frozen water?  FIXME
//
int set_pressure(const Epetra_MultiVector& pressure, double patm);


//
// Set the meteorological data.
// ------------------------------------------------------------------
// NOTE: all are of size ncolumns.
//
// NOTE: this is a little wonky, as CLM1D actually wants split precip,
//  and ATS uses split precip, but the CLM driver expects summed
//  precip that it splits by air temperature.  FIXME
//
// Input:
//   qSW        | Shortwave incoming radiation [W/m^2] Size ncolumns.
//   qLW        | Longwave incoming radiation [W/m^2] Size ncolumns.
//   pRain      | Rainfall precipitation rate [m/s]
//   pSnow      | Snowfall precipitation rate [m/s]
//   air_temp   | Air temperature [K]
//   vp_air     | Vapor pressure air [Pa]
//   wind_u     | Windspeed velocity [m/s]
//   p_atm      | Atmospheric pressure [Pa]
//
int set_met_data(const Epetra_MultiVector& qSW,
                 const Epetra_MultiVector& qLW,
                 const Epetra_MultiVector& pRain,
                 const Epetra_MultiVector& pSnow,
                 const Epetra_MultiVector& air_temp,
                 const Epetra_MultiVector& vp_air,
                 const Epetra_MultiVector& wind_u,
                 double p_atm);


//
// Advance the timestep
// ------------------------------------------------------------------
// Input:
//   step       | integer cycle number (logging only?)
//   time       | time at start of step (relative to zero time) [s]
//   dt         | step size [s]
int advance_time(int step, double time, double dt);


//
// Source/sink terms for integrated hydrology code.
// ------------------------------------------------------------------
// Output:
//   qW_surf    | water source/sink surface (sign?) [m/s] Size ncolumns.
//   qW_subsurf | water source/sink surface (sign?) [1/s] Size ncells.
//
int get_total_mass_fluxes(Epetra_MultiVector& qW_surf, Epetra_MultiVector& qW_subsurf);


//
// Energy fluxes for diagnostics/visualization
// ------------------------------------------------------------------
// Note: all sizes are ncolumns
//
// Output:
//   latent_heat        | Latent heat flux [W/m^2] (sign?)
//   sensible_heat      | Sensible heat flux [W/m^2] (sign?)
//   longwave_out       | Outward longwave radiation from surface to atmosphere [W/m^2]
//   conducted_e        | Energy conducted to subsurface (sign?) [W/m^2]
//
// The total variant includes canopy latent and sensible heat terms.  The
// ground variant does not.
int get_total_energy_fluxes(Epetra_MultiVector& latent_heat,
                            Epetra_MultiVector& sensible_heat,
                            Epetra_MultiVector& lw_out,
                            Epetra_MultiVector& conducted_e);
int get_ground_energy_fluxes(Epetra_MultiVector& latent_heat,
                             Epetra_MultiVector& sensible_heat,
                             Epetra_MultiVector& lw_out,
                             Epetra_MultiVector& conducted_e);

//
// Mass fluxes for mass balance and diagnostics/visualization
// ------------------------------------------------------------------
// Note: size ncolumns except where otherwise noted.
//
// Output:
//   evap_total         | Total evaporation to atmosphere [m/s]
//   evap_ground        | Ground component of evaporation (no snow sublimation) [m/s]
//   evap_soil          | Soil evaporation [m/s]
//   evap_canopy        | Canopy component of evaporation [m/s]
//   tran_veg           | Transpiration from vegetation over the column [m/s]
//   influx             | Precip? ???? [m/s]
//   irrigation         | Surface irrigation (can be intercepted) [m/s]
//   irrigation_inst    | Irrigation delivered to subsurface directly [1/s???] (Size ncells)
//   irrigation_flag    | ???
//   tran_soil          | Transpiration distributed to soil via rooting curve [1/s] (Size ncells)
//
int get_mass_fluxes(Epetra_MultiVector& evap_total,
                    Epetra_MultiVector& evap_ground,
                    Epetra_MultiVector& evap_soil,
                    Epetra_MultiVector& evap_canopy,
                    Epetra_MultiVector& tran_veg,
                    Epetra_MultiVector& influx,
                    Epetra_MultiVector& irrigation,
                    Epetra_MultiVector& inst_irrigation,
                    Epetra_MultiVector& irrigation_flag,
                    Epetra_MultiVector& tran_soil);


//
// Diagnostic variables
// ------------------------------------------------------------------
// Note: size ncolumns except where otherwise noted.
//
// Output:
//   swe                | snow water equivalent [m]
//   snow_depth         | [m] CHECK THIS!  CLM thinks this is [m] but everything
//                      |    else is [mm]... FIXME
//   canopy_storage     | water stored in the canopy [m]
//   T_skin             | Surface skin temperature [K]
//   T_veg              | Leaf temperature [K]
//   T_soil             | Soil temperature [K]  (size ncells)
//
int get_diagnostics(Epetra_MultiVector& swe,
                    Epetra_MultiVector& snow_depth,
                    Epetra_MultiVector& canopy_storage,
                    Epetra_MultiVector& Tskin,
                    Epetra_MultiVector& Tveg,
                    Epetra_MultiVector& Tsoil);


} // namespace CLM
} // namespace ATS


#endif
