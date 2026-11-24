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

#endif
