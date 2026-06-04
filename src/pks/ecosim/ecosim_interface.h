 /*****************************************************************************
 **
 ** C declarations of the ecosim interface
 **
 ******************************************************************************/

#include "data/BGC_containers.hh"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void ecosim_datatest();

void ecosim_setup(
  BGCProperties* properties,
  BGCState* state,
  BGCSizes* sizes,
  int num_iterations,
  int num_columns,
  int ncells_per_col_
);

void ecosim_shutdown();

void ecosim_advance(
  double delta_t,
  BGCProperties* properties,
  BGCState* state,
  BGCSizes* sizes,
  int num_iterations,
  int num_columns
);

#ifdef __cplusplus
}
#endif /* __cplusplus */
