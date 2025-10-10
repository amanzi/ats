/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
           Fenming Yuan
           Ethan Coon

This file defines a C interface to the ATS library
for use with LSMs.
*/
//! Wrapper for driving ATS from ELM.

#ifndef ELM_ATS_API_HH_
#define ELM_ATS_API_HH_

#ifdef __cplusplus
extern "C" {

// opaque pointer
// external caller only sees *ELM_ATSDriver_ptr - similar to void*, but better type safety
// ATS resolves ELM_ATSDriver_ptr as real ELM_ATSDriver during linking
class ELM_ATSDriver;
typedef ELM_ATSDriver *ELM_ATSDriver_ptr;

#else
// calling code should not dereference the pointer to the ATS object
// pointer hidden behind typedef to discourage
typedef struct ELM_ATSDriver_ptr *ELM_ATSDriver_ptr;

#endif

// allocate, call constructor and cast ptr to opaque ELM_ATSDriver_ptr
ELM_ATSDriver_ptr ats_create(MPI_Fint *f_comm, const char *input_filename);

// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATSDriver_ptr ats);


//
// Control
// -----------------------------------------------------------------------------

//
// These quantities should be compared against ELM to ensure consistent setup
//
void ats_get_mesh_info(ELM_ATSDriver_ptr ats,
                       int * const ncols_local,
                       int * const ncols_global,
                       int * const nlevgrnd,
                       double * const depth);

//
// simulation setup
//
void ats_setup(ELM_ATSDriver_ptr ats);

//
// initial conditions
//
void ats_initialize(ELM_ATSDriver_ptr ats);


// call driver advance(dt)
void ats_advance(ELM_ATSDriver_ptr ats,
                 double dt,
                 bool checkpoint,
                 bool visualize);

// call driver advance_test()
void ats_advance_test(ELM_ATSDriver_ptr ats);


//
// Memory movement/data passing
// -----------------------------------------------------------------------------
// double ats_get_scalar(ELM_ATSDriver_ptr ats, int scalar_id);
// void ats_set_scalar(ELM_ATSDriver_ptr ats, int scalar_id, double in);

void ats_get_field(ELM_ATSDriver_ptr ats, int var_id, double * const in);
double const * ats_get_field_ptr(ELM_ATSDriver_ptr ats, int var_id);
double * ats_get_field_ptr_w(ELM_ATSDriver_ptr ats, int var_id);
void ats_set_field(ELM_ATSDriver_ptr ats, int var_id, double const * const in);

#ifdef __cplusplus
}
#endif

// include guard
#endif
