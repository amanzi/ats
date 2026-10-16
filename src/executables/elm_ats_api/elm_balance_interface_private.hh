/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rich Fiorella, Ethan Coon
*/

//! Declaration of the ELM water-balance kernel bound from Fortran.
/*!

  elm_ats_water_balance_error_c is implemented in the shared, dependency-free
  Fortran module elm_water_balance_kernel.F90 (which lives in E3SM and is
  compiled into both the ELM build and this elm_ats_api library).  It returns
  ELM's column water-balance error errh2o [mm] given the subset of flux terms
  ATS controls.  All arguments are passed by reference (Fortran convention).

  Units (identical to ELM's BalanceCheckMod):
    endwb, begwb : column water storage end/begin of step [kg m-2 = mm]
    source       : net water source into the soil top             [mm s-1]
    evap, tran   : evaporation and transpiration sinks            [mm s-1]
    baseflow     : subsurface runoff (-> ELM qflx_drain)          [mm s-1]
    runoff       : surface runoff (-> ELM qflx_surf)              [mm s-1]
    dtime        : timestep                                        [s]
  Returns errh2o [mm].
*/

#ifndef ELM_BALANCE_INTERFACE_PRIVATE_HH_
#define ELM_BALANCE_INTERFACE_PRIVATE_HH_

#include <cstddef>

extern "C"
{
  double elm_ats_water_balance_error_c(const double* endwb,
                                       const double* begwb,
                                       const double* source,
                                       const double* evap,
                                       const double* tran,
                                       const double* baseflow,
                                       const double* runoff,
                                       const double* dtime);
}

inline double
elmWaterBalanceError(double endwb_mm, double begwb_mm, double source_mm_per_s,
                     double evap_mm_per_s, double tran_mm_per_s,
                     double baseflow_mm_per_s, double runoff_mm_per_s, double dt_s)
{
  return elm_ats_water_balance_error_c(&endwb_mm, &begwb_mm, &source_mm_per_s, &evap_mm_per_s,
                                       &tran_mm_per_s, &baseflow_mm_per_s, &runoff_mm_per_s,
                                       &dt_s);
}

//! Conversion from ATS storage [mol] to the ELM column storage [mm] the kernel wants.
/*!

  These are kept here, separate from ELM_ATSDriver, so the unit conversions can
  be tested without a mesh, State, or MPI -- the driver supplies Epetra row
  pointers and a mesh column view, while the tests supply plain arrays and a
  vector of indices.  Hence the templating: both are indexed with [] and expose
  size(), and nothing else is required of them.

  Molar density is per-cell because it is temperature/pressure dependent, so the
  division cannot be hoisted out of the sum.  The surface cell carries its own
  (surface) molar density, which is why it is a separate argument rather than
  another entry in the subsurface arrays.

*/

//! Water density [kg m-3].  ELM uses SHR_CONST_RHOFW from shr_const_mod.
constexpr double ELM_DENH2O = 1000.0;

//! m s-1 -> mm s-1, for the surface fluxes handed to the kernel.
constexpr double ELM_M_PER_S_TO_MM_PER_S = 1.0e3;

//! Total liquid water volume [m^3] of one column: its surface cell plus every
//! subsurface cell beneath it.  Summing both means the internal
//! surface<->subsurface exchange cancels, leaving only domain-boundary fluxes.
//!   surf_wc_mol [mol], surf_n_liq [mol m-3] -- the surface cell
//!   wc_mol [mol], n_liq [mol m-3]           -- indexed by subsurface cell id
//!   cells                                   -- this column's subsurface cell ids
template<class Vec, class Cells>
double
elmColumnWaterVolume(double surf_wc_mol, double surf_n_liq, const Vec& wc_mol, const Vec& n_liq,
                     const Cells& cells)
{
  double vol = surf_wc_mol / surf_n_liq;
  for (std::size_t j = 0; j != cells.size(); ++j) {
    const auto c = cells[j];
    vol += wc_mol[c] / n_liq[c];
  }
  return vol;
}

//! Column storage depth [kg m-2 = mm] from water volume [m^3] over the surface
//! cell area [m^2].  This is ELM's begwb/endwb.
inline double
elmStorageDepth(double vol_m3, double area_m2)
{
  return vol_m3 * ELM_DENH2O / area_m2;
}

#endif
