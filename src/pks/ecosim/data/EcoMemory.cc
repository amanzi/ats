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


/*******************************************************************************
 **
 **  Alquimia C memory utilities to handle memory management
 **
 **  Notes:
 **
 **  - calloc/malloc always return NULL pointers if they fail, so
 **    there is no need to pre-assign NULL for the pointers we are
 **    allocating here. For zero size or zero members, the returned
 **    pointer should be NULL or something that can be freed....
 **
 ** - free just releases the memory, it does not change the value
 **   of the pointer. After free, the pointer is no longer valid, so
 **   we set it to NULL.
 **
 *******************************************************************************/


#include "EcoMemory.hh"
#include "EcoContainers.hh"

// Returns the nearest power of 2 greater than or equal to n, or 0 if n == 0.
static inline int nearest_power_of_2(int n)
{
  if (n == 0) return 0;
  int twop = 1;
  while (twop < n)
    twop *= 2;
  return twop;
}

/*******************************************************************************
 **
 **  Vectors
 **
 *******************************************************************************/
void AllocateEcoVectorDouble(const int size, EcoVectorDouble* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (double*) calloc((size_t)vector->capacity, sizeof(double));
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}

void FreeEcoVectorDouble(EcoVectorDouble* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}

void AllocateEcoVectorInt(const int size, EcoVectorInt* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (int*) calloc((size_t)vector->capacity, sizeof(int));
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}

void FreeEcoVectorInt(EcoVectorInt* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}

/*******************************************************************************
 **
 **  Matrices
 **
 *******************************************************************************/
void AllocateEcoMatrixDouble(const int cells, const int columns, EcoMatrixDouble* matrix) {
  if ((cells > 0) || (columns > 0)) {
    matrix->cells = cells;
    matrix->columns = columns;
    matrix->capacity_cells = nearest_power_of_2(cells);
    matrix->capacity_columns = nearest_power_of_2(columns);
    matrix->data = (double*) calloc((size_t)matrix->capacity_cells * matrix->capacity_columns,
                                    sizeof(double));
  } else {
    matrix->cells = 0;
    matrix->columns = 0;
    matrix->capacity_cells = 0;
    matrix->capacity_columns = 0;
    matrix->data = NULL;
  }
}

void FreeEcoMatrixDouble(EcoMatrixDouble* matrix) {
  if (matrix != NULL) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->cells = 0;
    matrix->columns = 0;
    matrix->capacity_cells = 0;
    matrix->capacity_columns = 0;
  }
}

void AllocateEcoMatrixInt(const int cells, const int columns, EcoMatrixInt* matrix) {
  if ((cells > 0) || (columns > 0)) {
    matrix->cells = cells;
    matrix->columns = columns;
    matrix->capacity_cells = nearest_power_of_2(cells);
    matrix->capacity_columns = nearest_power_of_2(columns);
    matrix->data = (int*) calloc((size_t)matrix->capacity_cells * matrix->capacity_columns,
                                 sizeof(int));
  } else {
    matrix->cells = 0;
    matrix->columns = 0;
    matrix->capacity_cells = 0;
    matrix->capacity_columns = 0;
    matrix->data = NULL;
  }
}

void FreeEcoMatrixInt(EcoMatrixInt* matrix) {
  if (matrix != NULL) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->cells = 0;
    matrix->columns = 0;
    matrix->capacity_cells = 0;
    matrix->capacity_columns = 0;
  }
}

/*******************************************************************************
 **
 **  Tensors
 **
 *******************************************************************************/
void AllocateEcoTensorDouble(const int cells, const int columns, const int components,
                             EcoTensorDouble* tensor) {
  if ((cells > 0) || (columns > 0) || (components > 0)) {
    tensor->cells = cells;
    tensor->columns = columns;
    tensor->components = components;
    tensor->capacity_cells = nearest_power_of_2(cells);
    tensor->capacity_columns = nearest_power_of_2(columns);
    tensor->capacity_components = nearest_power_of_2(components);
    tensor->data = (double*) calloc((size_t)tensor->capacity_columns * tensor->capacity_cells *
                                      tensor->capacity_components,
                                    sizeof(double));
  } else {
    tensor->cells = 0;
    tensor->columns = 0;
    tensor->components = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
    tensor->data = NULL;
  }
}

void FreeEcoTensorDouble(EcoTensorDouble* tensor) {
  if (tensor != NULL) {
    free(tensor->data);
    tensor->data = NULL;
    tensor->cells = 0;
    tensor->columns = 0;
    tensor->components = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
  }
}

void AllocateEcoTensorInt(const int cells, const int columns, const int components,
                          EcoTensorInt* tensor) {
  if ((cells > 0) || (columns > 0) || (components > 0)) {
    tensor->cells = cells;
    tensor->columns = columns;
    tensor->components = components;
    tensor->capacity_cells = nearest_power_of_2(cells);
    tensor->capacity_columns = nearest_power_of_2(columns);
    tensor->capacity_components = nearest_power_of_2(components);
    tensor->data = (int*) calloc((size_t)tensor->capacity_columns * tensor->capacity_cells *
                                   tensor->capacity_components,
                                 sizeof(int));
  } else {
    tensor->cells = 0;
    tensor->columns = 0;
    tensor->components = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
    tensor->data = NULL;
  }
}

void FreeEcoTensorInt(EcoTensorInt* tensor) {
  if (tensor != NULL) {
    free(tensor->data);
    tensor->data = NULL;
    tensor->cells = 0;
    tensor->columns = 0;
    tensor->components = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
  }
}

/*******************************************************************************
 **
 **  Environment (ATS -> EcoSIM)
 **
 *******************************************************************************/
void AllocateEcoEnvironment(EcoEnvironment* env, int ncells_per_col_, int num_components,
                            int num_columns) {
  const int nc = ncells_per_col_, ncol = num_columns;
  AllocateEcoMatrixDouble(nc, ncol, &(env->liquid_density));
  AllocateEcoMatrixDouble(nc, ncol, &(env->gas_density));
  AllocateEcoMatrixDouble(nc, ncol, &(env->ice_density));
  AllocateEcoMatrixDouble(nc, ncol, &(env->rock_density));
  AllocateEcoMatrixDouble(nc, ncol, &(env->porosity));
  AllocateEcoMatrixDouble(nc, ncol, &(env->water_content));
  AllocateEcoMatrixDouble(nc, ncol, &(env->matric_pressure));
  AllocateEcoMatrixDouble(nc, ncol, &(env->temperature));
  AllocateEcoMatrixDouble(nc, ncol, &(env->hydraulic_conductivity));
  AllocateEcoMatrixDouble(nc, ncol, &(env->bulk_density));
  AllocateEcoMatrixDouble(nc, ncol, &(env->liquid_saturation));
  AllocateEcoMatrixDouble(nc, ncol, &(env->gas_saturation));
  AllocateEcoMatrixDouble(nc, ncol, &(env->ice_saturation));
  AllocateEcoMatrixDouble(nc, ncol, &(env->relative_permeability));
  AllocateEcoMatrixDouble(nc, ncol, &(env->thermal_conductivity));
  AllocateEcoMatrixDouble(nc, ncol, &(env->volume));
  AllocateEcoMatrixDouble(nc, ncol, &(env->depth));
  AllocateEcoMatrixDouble(nc, ncol, &(env->dz));
  AllocateEcoMatrixDouble(nc, ncol, &(env->plant_wilting_factor));
  AllocateEcoMatrixDouble(nc, ncol, &(env->rooting_depth_fraction));
  AllocateEcoMatrixDouble(nc, ncol, &(env->plant_functional_type));
  AllocateEcoMatrixDouble(nc, ncol, &(env->LAI));
  AllocateEcoMatrixDouble(nc, ncol, &(env->SAI));
  AllocateEcoTensorDouble(nc, ncol, num_components, &(env->mole_fraction));
  AllocateEcoVectorDouble(ncol, &(env->column_area));
  AllocateEcoVectorDouble(ncol, &(env->shortwave_radiation));
  AllocateEcoVectorDouble(ncol, &(env->longwave_radiation));
  AllocateEcoVectorDouble(ncol, &(env->air_temperature));
  AllocateEcoVectorDouble(ncol, &(env->vapor_pressure_air));
  AllocateEcoVectorDouble(ncol, &(env->wind_speed));
  AllocateEcoVectorDouble(ncol, &(env->precipitation));
  AllocateEcoVectorDouble(ncol, &(env->precipitation_snow));
  AllocateEcoVectorDouble(ncol, &(env->elevation));
  AllocateEcoVectorDouble(ncol, &(env->aspect));
  AllocateEcoVectorDouble(ncol, &(env->slope));
  AllocateEcoVectorDouble(ncol, &(env->vegetation_type));
  AllocateEcoVectorDouble(ncol, &(env->snow_albedo));
}

void FreeEcoEnvironment(EcoEnvironment* env) {
  if (env != NULL) {
    FreeEcoMatrixDouble(&(env->liquid_density));
    FreeEcoMatrixDouble(&(env->gas_density));
    FreeEcoMatrixDouble(&(env->ice_density));
    FreeEcoMatrixDouble(&(env->rock_density));
    FreeEcoMatrixDouble(&(env->porosity));
    FreeEcoMatrixDouble(&(env->water_content));
    FreeEcoMatrixDouble(&(env->matric_pressure));
    FreeEcoMatrixDouble(&(env->temperature));
    FreeEcoMatrixDouble(&(env->hydraulic_conductivity));
    FreeEcoMatrixDouble(&(env->bulk_density));
    FreeEcoMatrixDouble(&(env->liquid_saturation));
    FreeEcoMatrixDouble(&(env->gas_saturation));
    FreeEcoMatrixDouble(&(env->ice_saturation));
    FreeEcoMatrixDouble(&(env->relative_permeability));
    FreeEcoMatrixDouble(&(env->thermal_conductivity));
    FreeEcoMatrixDouble(&(env->volume));
    FreeEcoMatrixDouble(&(env->depth));
    FreeEcoMatrixDouble(&(env->dz));
    FreeEcoMatrixDouble(&(env->plant_wilting_factor));
    FreeEcoMatrixDouble(&(env->rooting_depth_fraction));
    FreeEcoMatrixDouble(&(env->plant_functional_type));
    FreeEcoTensorDouble(&(env->mole_fraction));
    FreeEcoVectorDouble(&(env->column_area));
    FreeEcoVectorDouble(&(env->shortwave_radiation));
    FreeEcoVectorDouble(&(env->longwave_radiation));
    FreeEcoVectorDouble(&(env->air_temperature));
    FreeEcoVectorDouble(&(env->vapor_pressure_air));
    FreeEcoVectorDouble(&(env->wind_speed));
    FreeEcoVectorDouble(&(env->precipitation));
    FreeEcoVectorDouble(&(env->precipitation_snow));
    FreeEcoVectorDouble(&(env->elevation));
    FreeEcoVectorDouble(&(env->aspect));
    FreeEcoVectorDouble(&(env->slope));
    FreeEcoMatrixDouble(&(env->LAI));
    FreeEcoMatrixDouble(&(env->SAI));
    FreeEcoVectorDouble(&(env->vegetation_type));
    FreeEcoVectorDouble(&(env->snow_albedo));
  }
}

/*******************************************************************************
 **
 **  Feedback (EcoSIM -> ATS)
 **
 *******************************************************************************/
void AllocateEcoFeedback(EcoFeedback* feedback, int ncells_per_col_, int num_columns, int num_pfts) {
  AllocateEcoMatrixDouble(ncells_per_col_, num_columns, &(feedback->subsurface_water_source));
  AllocateEcoMatrixDouble(ncells_per_col_, num_columns, &(feedback->subsurface_energy_source));
  AllocateEcoVectorDouble(num_columns, &(feedback->surface_water_source));
  AllocateEcoVectorDouble(num_columns, &(feedback->surface_energy_source));
  AllocateEcoVectorDouble(num_columns, &(feedback->snow_depth));
  AllocateEcoMatrixDouble(num_pfts, num_columns, &(feedback->canopy_snow));
}

void FreeEcoFeedback(EcoFeedback* feedback) {
  if (feedback != NULL) {
    FreeEcoMatrixDouble(&(feedback->subsurface_water_source));
    FreeEcoMatrixDouble(&(feedback->subsurface_energy_source));
    FreeEcoVectorDouble(&(feedback->surface_water_source));
    FreeEcoVectorDouble(&(feedback->surface_energy_source));
    FreeEcoVectorDouble(&(feedback->snow_depth));
    FreeEcoMatrixDouble(&(feedback->canopy_snow));
  }
}

/*******************************************************************************
 **
 **  EcoSIM internal state
 **
 *******************************************************************************/
void AllocateEcoInternalState(EcoInternalState* internal_state, int layout_version,
                              int num_entries, int num_columns, int values_per_column) {
  internal_state->layout_version = layout_version;
  internal_state->num_entries = num_entries;
  internal_state->num_columns = num_columns;
  internal_state->values_per_column = values_per_column;
  AllocateEcoMatrixDouble(values_per_column, num_columns, &(internal_state->values));
}

void FreeEcoInternalState(EcoInternalState* internal_state) {
  if (internal_state != NULL) {
    FreeEcoMatrixDouble(&(internal_state->values));
  }
}
