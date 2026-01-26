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

#include <iostream>
#include "BGC_memory.hh"
#include "BGC_containers.hh"

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
 **  BGC Vectors
 **
 *******************************************************************************/
void AllocateBGCVectorDouble(const int size, BGCVectorDouble* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (double*) calloc((size_t)vector->capacity, sizeof(double));
    //ALQUIMIA_ASSERT(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateBGCVectorDouble() */

void FreeBGCVectorDouble(BGCVectorDouble* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeBGCVectorDouble() */

void AllocateBGCVectorInt(const int size, BGCVectorInt* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (int*) calloc((size_t)vector->capacity, sizeof(int));
    //ALQUIMIA_ASSERT(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateBGCVectorInt() */

void FreeBGCVectorInt(BGCVectorInt* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeBGCVectorInt() */

void AllocateBGCVectorString(const int size, BGCVectorString* vector) {
  int i;
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (char**) calloc((size_t)vector->capacity, sizeof(char*));
    //ALQUIMIA_ASSERT(NULL != vector->data);
    for (i = 0; i < vector->size; ++i) {
      vector->data[i] = (char*) calloc((size_t)kBGCMaxStringLength, sizeof(char));
      //ALQUIMIA_ASSERT(NULL != vector->data[i]);
    }
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateBGCVectorString() */

void FreeBGCVectorString(BGCVectorString* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      free(vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeBGCVectorString() */

/*******************************************************************************
 **
 **  BGC Matrix
 **
 *******************************************************************************/
void AllocateBGCMatrixDouble(const int cells, const int columns, BGCMatrixDouble* matrix) {
  if ((cells > 0 ) || (columns > 0)){
    matrix->cells = cells;
    matrix->columns = columns;
    matrix->capacity_cells = nearest_power_of_2(cells);
    matrix->capacity_columns = nearest_power_of_2(columns);
    matrix->data = (double*) calloc((size_t)matrix->capacity_cells * matrix->capacity_columns, sizeof(double));

    //for (int i = 0; i < matrix->columns; ++i) {
    //  matrix->data[i] = (double*) calloc((size_t)matrix->capacity_cells, sizeof(double));
    //}
    //ALQUIMIA_ASSERT(NULL != matrix->data);
  } else {
    matrix->cells= 0;
    matrix->columns= 0;
    matrix->capacity_cells= 0;
    matrix->capacity_columns= 0;
    matrix->data = NULL;
  }
}  /* end AllocateBGCmatrixDouble() */

void FreeBGCMatrixDouble(BGCMatrixDouble* matrix) {
  if (matrix != NULL) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->cells= 0;
    matrix->columns= 0;
    matrix->capacity_cells= 0;
    matrix->capacity_columns= 0;
  }
}  /* end FreeBGCmatrixDouble() */

void AllocateBGCMatrixInt(const int cells, const int columns, BGCMatrixInt* matrix) {
  if ((cells> 0) || (columns> 0)) {
    matrix->cells= cells;
    matrix->columns= columns;
    matrix->capacity_cells= nearest_power_of_2(cells);
    matrix->capacity_columns= nearest_power_of_2(columns);
    matrix->data = (int*) calloc((size_t)matrix->capacity_cells * matrix->capacity_columns, sizeof(int));
    //for (int i = 0; i < matrix->columns; ++i) {
    //  matrix->data[i] = (int*) calloc((size_t)matrix->capacity_cells, sizeof(int));
    //}
    //ALQUIMIA_ASSERT(NULL != matrix->data);
  } else {
    matrix->cells= 0;
    matrix->columns= 0;
    matrix->capacity_cells= 0;
    matrix->capacity_columns= 0;
    matrix->data = NULL;
  }
}  /* end AllocateBGCMatrixInt() */

void FreeBGCMatrixInt(BGCMatrixInt* matrix) {
  if (matrix != NULL) {
    free(matrix->data);
    matrix->data = NULL;
    matrix->cells= 0;
    matrix->columns= 0;
    matrix->capacity_columns= 0;
    matrix->capacity_cells= 0;
  }
}  /* end FreeBGCMatrixInt() */

/*******************************************************************************
 **
 **  BGC Tensor
 **
 *******************************************************************************/

void AllocateBGCTensorDouble(const int cells, const int columns, const int components, BGCTensorDouble* tensor) {
  if ((cells> 0 ) || (columns> 0) || (components > 0)){
    tensor->cells = cells;
    tensor->columns = columns;
    tensor->components = components;

    tensor->capacity_cells= nearest_power_of_2(cells);
    tensor->capacity_columns= nearest_power_of_2(columns);
    tensor->capacity_components = nearest_power_of_2(components);

    tensor->data = (double***) calloc((size_t)tensor->capacity_columns, sizeof(double**));
    for (int i = 0; i < tensor->columns; ++i) {
      tensor->data[i] = (double**) calloc((size_t)tensor->capacity_cells, sizeof(double*));
      for (int j = 0; j < tensor->cells; ++j) {
        tensor->data[i][j] = (double*) calloc((size_t)tensor->capacity_components, sizeof(double));
      }
    }
    //ALQUIMIA_ASSERT(NULL != matrix->data);
  } else {
    tensor->cells= 0;
    tensor->columns= 0;
    tensor->components = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
    tensor->data = NULL;
  }
}  /* end AllocateBGCmatrixDouble() */

void FreeBGCTensorDouble(BGCTensorDouble* tensor) {
  if (tensor != NULL) {
    free(tensor->data);
    tensor->data = NULL;
    tensor->cells= 0;
    tensor->columns = 0;
    tensor->capacity_cells = 0;
    tensor->capacity_columns = 0;
    tensor->capacity_components = 0;
  }
}  /* end FreeBGCmatrixDouble() */

void AllocateBGCTensorInt(const int cells, const int columns, const int components, BGCTensorInt* tensor) {
  if ((cells> 0 ) || (columns> 0) || (components > 0)){
    tensor->cells= cells;
    tensor->columns= columns;
    tensor->components = components;

    tensor->capacity_cells= nearest_power_of_2(cells);
    tensor->capacity_columns= nearest_power_of_2(columns);
    tensor->capacity_components = nearest_power_of_2(components);

    tensor->data = (int***) calloc((size_t)tensor->capacity_columns, sizeof(int**));
    for (int i = 0; i < tensor->columns; ++i) {
      tensor->data[i] = (int**) calloc((size_t)tensor->capacity_cells, sizeof(int*));
      for (int j = 0; j < tensor->cells; ++j) {
        tensor->data[i][j] = (int*) calloc((size_t)tensor->capacity_components, sizeof(int));
      }
    }
    //ALQUIMIA_ASSERT(NULL != matrix->data);
  } else {
    tensor->cells= 0;
    tensor->columns= 0;
    tensor->components = 0;
    tensor->capacity_cells= 0;
    tensor->capacity_columns= 0;
    tensor->capacity_components = 0;
    tensor->data = NULL;
  }
}  /* end AllocateBGCmatrixint() */

void FreeBGCTensorInt(BGCTensorInt* tensor) {
  if (tensor != NULL) {
    free(tensor->data);
    tensor->data = NULL;
    tensor->cells= 0;
    tensor->columns= 0;
    tensor->capacity_cells= 0;
    tensor->capacity_columns= 0;
    tensor->capacity_components = 0;
  }
}  /* end FreeBGCmatrixint() */

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/
/*Note that sizes for all the datasets that are single vectors should just be
the size of the column, need to test
For reference the old function call was:
void AllocateBGCState(const BGCSizes* const sizes,
                           BGCState* state)*/

 void AllocateBGCState(BGCSizes* sizes, BGCState* state,
                       int ncells_per_col_, int num_components, int num_columns) {
   sizes->ncells_per_col_ = ncells_per_col_;
   sizes->num_components = num_components;
   sizes->num_columns = num_columns;

   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->liquid_density));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->gas_density));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->ice_density));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->rock_density));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->porosity));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->water_content));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->matric_pressure));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->temperature));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->hydraulic_conductivity));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->bulk_density));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->subsurface_energy_source));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(state->subsurface_water_source));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->surface_water_source));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->surface_energy_source));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->snow_depth));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->canopy_longwave_radiation));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->boundary_latent_heat_flux));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->boundary_sensible_heat_flux));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->canopy_surface_water));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->transpiration));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->evaporation_canopy));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->evaporation_bare_ground));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->evaporation_litter));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->evaporation_snow));
   AllocateBGCVectorDouble(sizes->num_columns, &(state->sublimation_snow));
   AllocateBGCMatrixDouble(sizes->num_columns, sizes->num_columns, &(state->snow_temperature));
   AllocateBGCTensorDouble(sizes->ncells_per_col_, sizes->num_columns, sizes->num_components, &(state->total_component_concentration));
   //ALQUIMIA_ASSERT(state->total_mobile.data != NULL);
 }  /* end AllocateBGCState() */

 void FreeBGCState(BGCState* state) {
   if (state != NULL) {
     FreeBGCMatrixDouble(&(state->liquid_density));
     FreeBGCMatrixDouble(&(state->gas_density));
     FreeBGCMatrixDouble(&(state->ice_density));
     FreeBGCMatrixDouble(&(state->rock_density));
     FreeBGCMatrixDouble(&(state->porosity));
     FreeBGCMatrixDouble(&(state->water_content));
     FreeBGCMatrixDouble(&(state->matric_pressure));
     FreeBGCMatrixDouble(&(state->temperature));
     FreeBGCMatrixDouble(&(state->hydraulic_conductivity));
     FreeBGCMatrixDouble(&(state->bulk_density));
     FreeBGCMatrixDouble(&(state->subsurface_energy_source));
     FreeBGCMatrixDouble(&(state->subsurface_water_source));
     FreeBGCVectorDouble(&(state->surface_energy_source));
     FreeBGCVectorDouble(&(state->surface_water_source));
     FreeBGCVectorDouble(&(state->snow_depth));
     FreeBGCVectorDouble(&(state->canopy_longwave_radiation));
     FreeBGCVectorDouble(&(state->boundary_latent_heat_flux));
     FreeBGCVectorDouble(&(state->boundary_sensible_heat_flux));
     FreeBGCVectorDouble(&(state->canopy_surface_water));
     FreeBGCVectorDouble(&(state->transpiration));
     FreeBGCVectorDouble(&(state->evaporation_canopy));
     FreeBGCVectorDouble(&(state->evaporation_bare_ground));
     FreeBGCVectorDouble(&(state->evaporation_litter));
     FreeBGCVectorDouble(&(state->evaporation_snow));
     FreeBGCVectorDouble(&(state->sublimation_snow));
     FreeBGCMatrixDouble(&(state->snow_temperature));
     FreeBGCTensorDouble(&(state->total_component_concentration));
   }
 }  /* end FreeAlquimiaState() */

 /*******************************************************************************
  **
  **  Auxiliary Data
  **
  *******************************************************************************/
 /*
 void AllocateBGCAuxiliaryData(const BGCSizes* const sizes, BGCAuxiliaryData* aux_data,
                               int ncells_per_col_) {
   AllocateBGCMatrixInt(sizes->ncells_per_col_,
                             &(aux_data->aux_ints));

   AllocateBGCMatrixDouble(sizes->ncells_per_col_,
                                &(aux_data->aux_doubles));

 }  // end AllocateAlquimiaAuxiliaryData()

 void FreeBGCAuxiliaryData(BGCAuxiliaryData* aux_data) {
   if (aux_data != NULL) {
     FreeBGCMatrixInt(&(aux_data->aux_ints));
     FreeBGCMatrixDouble(&(aux_data->aux_doubles));
   }
 }  // end FreeAlquimiaAuxiliaryData()
 */

 /*******************************************************************************
  **
  **  Properties
  **
  *******************************************************************************/

 void AllocateBGCProperties(BGCSizes* sizes, BGCProperties* properties,
                           int ncells_per_col_, int num_columns) {

   sizes->ncells_per_col_ = ncells_per_col_;
   sizes->num_columns = num_columns;
   //sizes->num_components = num_components;
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->liquid_saturation));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->gas_saturation));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->ice_saturation));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->relative_permeability));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->thermal_conductivity));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->volume));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->depth));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->depth_c));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->dz));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->plant_wilting_factor));
   AllocateBGCMatrixDouble(sizes->ncells_per_col_, sizes->num_columns, &(properties->rooting_depth_fraction));

   AllocateBGCVectorDouble(sizes->num_columns, &(properties->column_area));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->shortwave_radiation));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->longwave_radiation));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->air_temperature));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->vapor_pressure_air));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->wind_speed));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->precipitation));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->precipitation_snow));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->elevation));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->aspect));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->slope));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->LAI));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->SAI));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->vegetation_type));
   AllocateBGCVectorDouble(sizes->num_columns, &(properties->snow_albedo));
 }  /* end AllocateAlquimiaProperties() */

 void FreeBGCProperties(BGCProperties* properties) {
   if (properties != NULL) {
     FreeBGCMatrixDouble(&(properties->liquid_saturation));
     FreeBGCMatrixDouble(&(properties->gas_saturation));
     FreeBGCMatrixDouble(&(properties->ice_saturation));
     FreeBGCMatrixDouble(&(properties->relative_permeability));
     FreeBGCMatrixDouble(&(properties->thermal_conductivity));
     FreeBGCMatrixDouble(&(properties->volume));
     FreeBGCMatrixDouble(&(properties->depth));
     FreeBGCMatrixDouble(&(properties->dz));
     FreeBGCMatrixDouble(&(properties->plant_wilting_factor));
     FreeBGCMatrixDouble(&(properties->rooting_depth_fraction));

     FreeBGCVectorDouble(&(properties->column_area));
     FreeBGCVectorDouble(&(properties->shortwave_radiation));
     FreeBGCVectorDouble(&(properties->longwave_radiation));
     FreeBGCVectorDouble(&(properties->air_temperature));
     FreeBGCVectorDouble(&(properties->vapor_pressure_air));
     FreeBGCVectorDouble(&(properties->wind_speed));
     FreeBGCVectorDouble(&(properties->precipitation));
     FreeBGCVectorDouble(&(properties->precipitation_snow));
     FreeBGCVectorDouble(&(properties->elevation));
     FreeBGCVectorDouble(&(properties->aspect));
     FreeBGCVectorDouble(&(properties->slope));
     FreeBGCVectorDouble(&(properties->LAI));
     FreeBGCVectorDouble(&(properties->SAI));
     FreeBGCVectorDouble(&(properties->vegetation_type));
     FreeBGCVectorDouble(&(properties->snow_albedo));
   }
 }

/* OLD VERSION OF THE DATA MODULES
void AllocateBGCState(BGCSizes* sizes, BGCState* state,
                      int ncells_per_col_, int num_components) {
  sizes->ncells_per_col_ = ncells_per_col_;
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->liquid_density));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->gas_density));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->ice_density));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->porosity));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->water_content));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->suction_head));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->temperature));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->hydraulic_conductivity));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(state->bulk_density));
  AllocateBGCMatrixDouble(sizes->ncells_per_col_,sizes->num_components, &(state->total_component_concentration));
  //ALQUIMIA_ASSERT(state->total_mobile.data != NULL);

}

void FreeBGCState(BGCState* state) {
  if (state != NULL) {
    FreeBGCVectorDouble(&(state->liquid_density));
    FreeBGCVectorDouble(&(state->gas_density));
    FreeBGCVectorDouble(&(state->ice_density));
    FreeBGCVectorDouble(&(state->porosity));
    FreeBGCVectorDouble(&(state->water_content));
    FreeBGCVectorDouble(&(state->suction_head));
    FreeBGCVectorDouble(&(state->temperature));
    FreeBGCVectorDouble(&(state->hydraulic_conductivity));
    FreeBGCVectorDouble(&(state->bulk_density));
    FreeBGCMatrixDouble(&(state->total_component_concentration));
  }
}

void AllocateBGCAuxiliaryData(const BGCSizes* const sizes, BGCAuxiliaryData* aux_data,
                              int ncells_per_col_) {
  AllocateBGCVectorInt(sizes->ncells_per_col_,
                            &(aux_data->aux_ints));

  AllocateBGCVectorDouble(sizes->ncells_per_col_,
                               &(aux_data->aux_doubles));

}

void FreeBGCAuxiliaryData(BGCAuxiliaryData* aux_data) {
  if (aux_data != NULL) {
    FreeBGCVectorInt(&(aux_data->aux_ints));
    FreeBGCVectorDouble(&(aux_data->aux_doubles));
  }
}

void AllocateBGCProperties(BGCSizes* sizes, BGCProperties* properties,
                          int ncells_per_col_) {
  sizes->ncells_per_col_ = ncells_per_col_;

  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->liquid_saturation));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->gas_saturation));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->ice_saturation));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->relative_permeability));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->thermal_conductivity));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->volume));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->depth));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->dz));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->plant_wilting_factor));
  AllocateBGCVectorDouble(sizes->ncells_per_col_, &(properties->rooting_depth_fraction));
}

void FreeBGCProperties(BGCProperties* properties) {
  if (properties != NULL) {
    FreeBGCVectorDouble(&(properties->liquid_saturation));
    FreeBGCVectorDouble(&(properties->gas_saturation));
    FreeBGCVectorDouble(&(properties->ice_saturation));
    FreeBGCVectorDouble(&(properties->relative_permeability));
    FreeBGCVectorDouble(&(properties->thermal_conductivity));
    FreeBGCVectorDouble(&(properties->volume));
    FreeBGCVectorDouble(&(properties->depth));
    FreeBGCVectorDouble(&(properties->dz));
    FreeBGCVectorDouble(&(properties->plant_wilting_factor));
    FreeBGCVectorDouble(&(properties->rooting_depth_fraction));
  }
}

/*******************************************************************************
 **
 **  Problem Meta Data
 **
 *******************************************************************************/
/*
void AllocateAlquimiaProblemMetaData(const AlquimiaSizes* const sizes,
                                     AlquimiaProblemMetaData* meta_data) {

  AllocateAlquimiaVectorString(sizes->num_primary, &(meta_data->primary_names));
  ALQUIMIA_ASSERT(meta_data->primary_names.data != NULL);

  AllocateAlquimiaVectorInt(sizes->num_primary, &(meta_data->positivity));
  memset(meta_data->positivity.data, 0, sizeof(int) * sizes->num_primary);

  AllocateAlquimiaVectorString(sizes->num_minerals,
                               &(meta_data->mineral_names));

  AllocateAlquimiaVectorString(sizes->num_surface_sites,
                               &(meta_data->surface_site_names));

  AllocateAlquimiaVectorString(sizes->num_ion_exchange_sites,
                               &(meta_data->ion_exchange_names));

  AllocateAlquimiaVectorString(sizes->num_isotherm_species,
                               &(meta_data->isotherm_species_names));

  AllocateAlquimiaVectorString(sizes->num_aqueous_kinetics,
                               &(meta_data->aqueous_kinetic_names));

}  //end AllocateAlquimiaProblemMetaData()

void FreeAlquimiaProblemMetaData(AlquimiaProblemMetaData* meta_data) {

  if (meta_data != NULL) {
    FreeAlquimiaVectorString(&(meta_data->primary_names));
    FreeAlquimiaVectorInt(&(meta_data->positivity));
    FreeAlquimiaVectorString(&(meta_data->mineral_names));
    FreeAlquimiaVectorString(&(meta_data->surface_site_names));
    FreeAlquimiaVectorString(&(meta_data->ion_exchange_names));
    FreeAlquimiaVectorString(&(meta_data->isotherm_species_names));
    FreeAlquimiaVectorString(&(meta_data->aqueous_kinetic_names));
  }
}  end FreeAlquimiaProblemMetaData() */

/*******************************************************************************
 **
 **  Auxiliary Output Data
 **
 *******************************************************************************/
/*
void AllocateAlquimiaAuxiliaryOutputData(const AlquimiaSizes* const sizes,
                                         AlquimiaAuxiliaryOutputData* aux_output) {
  aux_output->pH = -999.9;
  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(aux_output->mineral_saturation_index));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_kinetics,
                               &(aux_output->aqueous_kinetic_rate));

  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(aux_output->mineral_reaction_rate));

  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_activity_coeff));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_activity_coeff));

}  // end AllocateAlquimiaAuxiliaryOutputData()

void FreeAlquimiaAuxiliaryOutputData(AlquimiaAuxiliaryOutputData* aux_output) {
  if (aux_output != NULL) {
    FreeAlquimiaVectorDouble(&(aux_output->aqueous_kinetic_rate));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_saturation_index));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_reaction_rate));
    FreeAlquimiaVectorDouble(&(aux_output->primary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->primary_activity_coeff));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_activity_coeff));
  }
}   end FreeAlquimiaAuxiliaryOutputData() */

/*******************************************************************************
 **
 **  Engine Status
 **
 *******************************************************************************/
/*
void AllocateAlquimiaEngineStatus(AlquimiaEngineStatus* status) {

  status->message = (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  if (NULL == status->message) {
    // TODO(bja): error handling
  }
}  // end AllocateAlquimiaEngineStatus()

void FreeAlquimiaEngineStatus(AlquimiaEngineStatus* status) {
  if (status != NULL) {
    free(status->message);
  }
  status->message = NULL;

}  end FreeAlquimiaEngineStatus() */



/*******************************************************************************
 **
 **  Geochemical conditions/constraints
 **
 *******************************************************************************/
/*
void AllocateAlquimiaGeochemicalConditionVector(const int num_conditions,
                                                AlquimiaGeochemicalConditionVector* condition_list) {
  // NOTE: we are only allocating pointers to N conditions here, not
     the actual conditions themselves.
  fprintf(stdout, " AllocateAlquimiaGeochemicalConditionList() : %d\n",
          num_conditions);
  condition_list->size = num_conditions;
  condition_list->capacity = nearest_power_of_2(num_conditions);

  if (condition_list->size > 0) {
    condition_list->data = (AlquimiaGeochemicalCondition*)
        calloc((size_t)condition_list->capacity,
               sizeof(AlquimiaGeochemicalCondition));
  }
}  // end AllocateAlquimiaGeochemicalConditionVector()

void AllocateAlquimiaGeochemicalCondition(const int size_name,
                                          const int num_aqueous_constraints,
                                          const int num_mineral_constraints,
    AlquimiaGeochemicalCondition* condition) {
  // NOTE: we are only allocating pointers to N constraints here, not
     the actual condstraints themselves.
  if (condition != NULL) {
    // size_name + 1 to include the null character!
    condition->name = (char*) calloc((size_t)size_name+1, sizeof(char));
    AllocateAlquimiaAqueousConstraintVector(num_aqueous_constraints, &condition->aqueous_constraints);
    AllocateAlquimiaMineralConstraintVector(num_mineral_constraints, &condition->mineral_constraints);
  }
}  // end AllocateAlquimiaGeochemicalCondition()

void AllocateAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint) {
  constraint->primary_species_name =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->constraint_type =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->associated_species =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->value = 0.0;
}  // end AllocateAlquimiaAqueousConstraint()

void AllocateAlquimiaAqueousConstraintVector(int num_constraints,
                                             AlquimiaAqueousConstraintVector* constraint_list) {
  constraint_list->size = num_constraints;
  constraint_list->capacity = nearest_power_of_2(num_constraints);
  if (constraint_list->size > 0) {
    constraint_list->data = (AlquimiaAqueousConstraint*)
      calloc((size_t)constraint_list->capacity,
             sizeof(AlquimiaAqueousConstraint));
  }
  else
    constraint_list->data = NULL;
}

void AllocateAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint) {
  constraint->mineral_name =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->volume_fraction = -1.0;
  constraint->specific_surface_area = -1.0;
}  // end AllocateAlquimiaMineralConstraint()

void AllocateAlquimiaMineralConstraintVector(int num_constraints,
                                             AlquimiaMineralConstraintVector* constraint_list) {
  constraint_list->size = num_constraints;
  constraint_list->capacity = nearest_power_of_2(num_constraints);
  if (constraint_list->size > 0) {
    constraint_list->data = (AlquimiaMineralConstraint*)
      calloc((size_t)constraint_list->capacity,
             sizeof(AlquimiaMineralConstraint));
  }
  else
    constraint_list->data = NULL;
}

void FreeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* condition_list) {
  int i;
  if (condition_list != NULL) {
    for (i = 0; i < condition_list->size; ++i) {
      FreeAlquimiaGeochemicalCondition(&(condition_list->data[i]));
    }
    if (condition_list->data != NULL) {
      free(condition_list->data);
      condition_list->data = NULL;
    }
    condition_list->size = 0;
    condition_list->capacity = 0;
  }
}  // end FreeAlquimiaGeochemicalConditionList()

void FreeAlquimiaGeochemicalCondition(AlquimiaGeochemicalCondition* condition) {
  if (condition != NULL) {
    if (condition->name != NULL) {
      free(condition->name);
      condition->name = NULL;
    }
    FreeAlquimiaAqueousConstraintVector(&(condition->aqueous_constraints));
    FreeAlquimiaMineralConstraintVector(&(condition->mineral_constraints));
  }
}  // end FreeAlquimiaGeochemicalCondition()

void FreeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaAqueousConstraint(&vector->data[i]);
    }
    if (vector->data != NULL) {
      free(vector->data);
      vector->data = NULL;
    }
    vector->size = 0;
    vector->capacity = 0;
  }
}  // end FreeAlquimiaAqueousConstraintVector()

void FreeAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint) {
  free(constraint->primary_species_name);
  constraint->primary_species_name = NULL;
  free(constraint->constraint_type);
  constraint->constraint_type = NULL;
  free(constraint->associated_species);
  constraint->associated_species = NULL;
}  // end FreeAlquimiaAqueousConstraint()

void FreeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaMineralConstraint(&vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  // end FreeAlquimiaMineralConstraintVector()

void FreeAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint) {
  free(constraint->mineral_name);
  constraint->mineral_name = NULL;
}  // end FreeAlquimiaMineralConstraint() */


/*******************************************************************************
 **
 **  Data convenience struct
 **
 *******************************************************************************/
/*
void AllocateAlquimiaData(AlquimiaData* data) {
    AllocateAlquimiaState(&data->sizes, &data->state);
    AllocateAlquimiaProperties(&data->sizes, &data->properties);
    AllocateAlquimiaAuxiliaryData(&data->sizes, &data->aux_data);
    AllocateAlquimiaProblemMetaData(&data->sizes, &data->meta_data);
    AllocateAlquimiaAuxiliaryOutputData(&data->sizes, &data->aux_output);
}  // end AllocateAlquimiaData()


void FreeAlquimiaData(AlquimiaData* data) {
  FreeAlquimiaState(&data->state);
  FreeAlquimiaProperties(&data->properties);
  FreeAlquimiaAuxiliaryData(&data->aux_data);
  FreeAlquimiaProblemMetaData(&data->meta_data);
  FreeAlquimiaAuxiliaryOutputData(&data->aux_output);
}  // end FreeAlquimiaData() */
