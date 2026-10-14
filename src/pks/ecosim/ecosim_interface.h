 /*****************************************************************************
 **
 ** C declarations of the ecosim interface (implemented in EcoSIM's
 ** ATSUtils/ecosim_wrappers.F90)
 **
 ******************************************************************************/

#ifndef ECOSIM_INTERFACE_H_
#define ECOSIM_INTERFACE_H_

#include "data/EcoContainers.hh"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void ecosim_datatest();

/* config is used only here; environment and feedback (snow depth) are also
   read here, before EcoSIM initializes */
void ecosim_setup(
  EcoConfig* config,
  EcoEnvironment* environment,
  EcoFeedback* feedback,
  EcoInternalState* internal_state,
  EcoSizes* sizes,
  int num_iterations,
  int num_columns,
  int ncells_per_col_
);

void ecosim_shutdown();

void ecosim_advance(
  double delta_t,
  EcoEnvironment* environment,
  EcoFeedback* feedback,
  EcoInternalState* internal_state,
  EcoSizes* sizes,
  int num_iterations,
  int num_columns
);

/* sizes in bytes of the bind(C) container types as EcoSIM sees them, in the
   order: EcoVectorDouble, EcoMatrixDouble, EcoTensorDouble, EcoSizes,
   EcoConfig, EcoEnvironment, EcoFeedback, EcoInternalState.
   n is kEcoNumContainerTypes on input, the number EcoSIM knows on output. */
void ecosim_container_sizes(int* n, size_t* sizes);

/* layout of the EcoSIM internal state; static, callable before ecosim_setup */
int ecosim_internal_state_layout_version();

int ecosim_internal_state_num_entries(const EcoSizes* sizes);

/* i is 0-based; buffers must hold kEcoStateNameLength, kEcoStateNameLength,
   kEcoStateUnitsLength and kEcoStateDescriptionLength characters */
void ecosim_internal_state_entry(int i,
                                 const EcoSizes* sizes,
                                 char* ats_name,
                                 char* ecosim_name,
                                 char* units,
                                 char* description,
                                 int* num_components,
                                 int* role);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ECOSIM_INTERFACE_H_ */
