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


#ifndef ECO_C_MEMORY_H_
#define ECO_C_MEMORY_H_

#include "EcoContainers.hh"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* Vectors */
  void AllocateEcoVectorDouble(const int size, EcoVectorDouble* vector);
  void FreeEcoVectorDouble(EcoVectorDouble* vector);

  void AllocateEcoVectorInt(const int size, EcoVectorInt* vector);
  void FreeEcoVectorInt(EcoVectorInt* vector);

  /* Matrices */
  void AllocateEcoMatrixDouble(const int cells, const int columns, EcoMatrixDouble* matrix);
  void FreeEcoMatrixDouble(EcoMatrixDouble* matrix);

  void AllocateEcoMatrixInt(const int cells, const int columns, EcoMatrixInt* matrix);
  void FreeEcoMatrixInt(EcoMatrixInt* matrix);

  /* Tensors */
  void AllocateEcoTensorDouble(const int cells, const int columns, const int components,
                               EcoTensorDouble* tensor);
  void FreeEcoTensorDouble(EcoTensorDouble* tensor);

  void AllocateEcoTensorInt(const int cells, const int columns, const int components,
                            EcoTensorInt* tensor);
  void FreeEcoTensorInt(EcoTensorInt* tensor);

  /* Exchange containers (EcoConfig holds no arrays and needs no allocation) */
  void AllocateEcoEnvironment(EcoEnvironment* environment,
                              int ncells_per_col_,
                              int num_components,
                              int num_columns);
  void FreeEcoEnvironment(EcoEnvironment* environment);

  void AllocateEcoFeedback(EcoFeedback* feedback, int ncells_per_col_, int num_columns, int num_pfts);
  void FreeEcoFeedback(EcoFeedback* feedback);

  void AllocateEcoInternalState(EcoInternalState* internal_state,
                                int layout_version,
                                int num_entries,
                                int num_columns,
                                int values_per_column);
  void FreeEcoInternalState(EcoInternalState* internal_state);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ECO_C_MEMORY_H_ */
