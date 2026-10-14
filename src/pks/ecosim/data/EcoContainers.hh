/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
** through Lawrence Berkeley National Laboratory (subject to receipt of any
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
**
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software,
** please contact Berkeley Lab's Technology Transfer and Intellectual Property
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
**
** NOTICE.  This software was developed under funding from the U.S. Department
** of Energy.  As such, the U.S. Government has been granted for itself and
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
** license in the Software to reproduce, prepare derivative works, and perform
** publicly and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of Energy,
** and subject to any subsequent five (5) year renewals, the U.S. Government is
** granted for itself and others acting on its behalf a paid-up, nonexclusive,
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly,
** and to permit others to do so.
**
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#ifndef ECO_CONTAINERS_H_
#define ECO_CONTAINERS_H_

/*******************************************************************************
 **
 ** Containers exchanged between ATS and EcoSIM (adapted from Alquimia).
 **
 ** These are passed directly to the Fortran routines, which mirror them as
 ** bind(C) types in EcoSIM's EcoContainers.F90 (EcoContainers_module). Member
 ** order and types must match exactly; ecosim_container_sizes() lets ATS
 ** check the sizes at setup.
 **
 ** Roles, by direction and lifetime:
 **   EcoSizes          dimensions
 **   EcoConfig         ATS -> EcoSIM, once at setup: run parameters and flags
 **   EcoEnvironment    ATS -> EcoSIM, at setup and every advance: soil state
 **                     and properties, geometry, forcing, prescribed
 **                     vegetation, clock. Owned by ATS; never read back.
 **   EcoFeedback       EcoSIM -> ATS, every advance: what ATS physics uses
 **                     from EcoSIM. snow_depth is also sent in (EcoSIM's snow
 **                     state); the incoming sources are not used by EcoSIM.
 **   EcoInternalState  EcoSIM-private carried state and outputs, stored by
 **                     ATS without interpreting it (ATSStateRegistryMod).
 **
 ** Members marked "not filled" or "placeholder" are carried over unchanged
 ** from the earlier BGCState/BGCProperties and are to be fixed separately.
 **
 ******************************************************************************/

#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct {
    int size, capacity;
    double* data;
  } EcoVectorDouble;

  typedef struct {
    int size, capacity;
    int* data;
  } EcoVectorInt;

  /* cells x columns, column-major: data[column*cells + cell] */
  typedef struct {
    int cells, columns, capacity_cells, capacity_columns;
    double* data;
  } EcoMatrixDouble;

  typedef struct {
    int cells, columns, capacity_cells, capacity_columns;
    int* data;
  } EcoMatrixInt;

  /* cells x columns x components */
  typedef struct {
    int cells, columns, components, capacity_cells, capacity_columns, capacity_components;
    double* data;
  } EcoTensorDouble;

  typedef struct {
    int cells, columns, components, capacity_cells, capacity_columns, capacity_components;
    int* data;
  } EcoTensorInt;

  typedef struct {
    int ncells_per_col_;
    int num_components;
    int num_columns;
    int num_pfts;
  } EcoSizes;

  /* ATS -> EcoSIM, once at setup */
  typedef struct {
    double heat_capacity;   /* [MJ mol^-1 K^-1] */
    double field_capacity;  /* pressure at field capacity [MPa] */
    double wilting_point;   /* pressure at wilting point [MPa] */
    bool p_bool;            /* EcoSIM precipitation (total, partitioned by EcoSIM) */
    bool a_bool;            /* prescribe snow albedo */
    bool pheno_bool;        /* prescribe phenology */
    bool microbe_bool;      /* microbe model (not copied by EcoSIM) */
    char* pft_file;         /* PFT parameter file, owned by the PK */
  } EcoConfig;

  /* ATS -> EcoSIM, at setup and every advance; never read back */
  typedef struct {
    /* per cell: ncells_per_col_ x num_columns */
    EcoMatrixDouble liquid_density;         /* molar_density_liquid */
    EcoMatrixDouble gas_density;            /* not filled */
    EcoMatrixDouble ice_density;            /* zeros when has_ice (not filled from a field) */
    EcoMatrixDouble rock_density;
    EcoMatrixDouble porosity;
    EcoMatrixDouble water_content;
    EcoMatrixDouble matric_pressure;        /* from capillary_pressure_gas_liq */
    EcoMatrixDouble temperature;
    EcoMatrixDouble hydraulic_conductivity; /* not used by EcoSIM */
    EcoMatrixDouble bulk_density;           /* not filled */
    EcoMatrixDouble liquid_saturation;
    EcoMatrixDouble gas_saturation;         /* not filled */
    EcoMatrixDouble ice_saturation;         /* zeros when has_ice (not filled from a field) */
    EcoMatrixDouble relative_permeability;  /* zeros (not filled from a field) */
    EcoMatrixDouble thermal_conductivity;   /* not filled */
    EcoMatrixDouble volume;
    EcoMatrixDouble depth;
    EcoMatrixDouble dz;
    EcoMatrixDouble plant_wilting_factor;   /* placeholder, filled from porosity */
    EcoMatrixDouble rooting_depth_fraction; /* placeholder, filled from porosity */
    EcoMatrixDouble plant_functional_type;  /* first num_pfts entries of each column */
    EcoMatrixDouble LAI;  /* num_pfts x num_columns */
    EcoMatrixDouble SAI;  /* num_pfts x num_columns */
    EcoTensorDouble mole_fraction;          /* ncells x num_columns x num_components; microbes only */
    /* per column: num_columns */
    EcoVectorDouble column_area;
    EcoVectorDouble shortwave_radiation;
    EcoVectorDouble longwave_radiation;     /* not filled */
    EcoVectorDouble air_temperature;
    EcoVectorDouble vapor_pressure_air;
    EcoVectorDouble wind_speed;
    EcoVectorDouble precipitation;          /* total if p_bool, else rain */
    EcoVectorDouble precipitation_snow;     /* only if !p_bool */
    EcoVectorDouble elevation;
    EcoVectorDouble aspect;
    EcoVectorDouble slope;
    EcoVectorDouble vegetation_type;        /* not filled */
    EcoVectorDouble snow_albedo;
    /* atmosphere composition (not initialized by the PK) */
    double atm_n2;
    double atm_o2;
    double atm_co2;
    double atm_ch4;
    double atm_n2o;
    double atm_h2;
    double atm_nh3;
    /* clock (set before every advance) */
    int current_day;
    int current_year;
  } EcoEnvironment;

  /* EcoSIM -> ATS, every advance */
  typedef struct {
    EcoMatrixDouble subsurface_water_source;  /* EcoSIM: m3 per grid cell per hour */
    EcoMatrixDouble subsurface_energy_source; /* not read back by ATS */
    EcoVectorDouble surface_water_source;     /* EcoSIM: m h-1 */
    EcoVectorDouble surface_energy_source;    /* EcoSIM: per hour */
    EcoVectorDouble snow_depth;               /* also sent in: EcoSIM's snow state */
    EcoMatrixDouble canopy_snow;              /* num_pfts x num_columns; also sent in */
  } EcoFeedback;

  /* EcoSIM-private data (carried state and EcoSIM-only outputs), packed by
     EcoSIM's ATSStateRegistryMod. ATS stores it without interpreting it; the
     layout is queried with ecosim_internal_state_entry(). */
  typedef struct {
    int layout_version;
    int num_entries;
    int num_columns;
    int values_per_column;
    EcoMatrixDouble values; /* values_per_column x num_columns */
  } EcoInternalState;

  /* roles of internal state entries, and string lengths of their layout */
  enum { kEcoRolePrivate = 0, kEcoRoleOutput = 1 };
  enum { kEcoStateNameLength = 64, kEcoStateUnitsLength = 32, kEcoStateDescriptionLength = 128 };

  /* container types whose sizes are checked against EcoSIM (same order as
     ecosim_container_sizes() in ecosim_wrappers.F90) */
  enum { kEcoNumContainerTypes = 8 };

  typedef struct {
    void (*DataTest)();

    void (*Setup)(
      EcoConfig* config,
      EcoEnvironment* environment,
      EcoFeedback* feedback,
      EcoInternalState* internal_state,
      EcoSizes* sizes,
      int num_iterations,
      int num_columns,
      int ncells_per_col_);

    void (*Shutdown)();

    void (*Advance)(
      double delta_t,
      EcoEnvironment* environment,
      EcoFeedback* feedback,
      EcoInternalState* internal_state,
      EcoSizes* sizes,
      int num_iterations,
      int num_columns);

  } EcoInterface;

  void CreateEcoInterface(const char* const engine_name, EcoInterface* interface);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ECO_CONTAINERS_H_ */
