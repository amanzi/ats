/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
           Fenming Yuan
           Ethan Coon
*/
//! Wrapper for driving ATS from ELM.

#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

#ifdef __cplusplus

namespace {
  // did ATS initialize kokkos?
  bool ats_kokkos_init = false;
}

extern "C"
{

// allocate, call constructor and cast ptr to opaque ELM_ATSDriver_ptr
ELM_ATSDriver_ptr ats_create(MPI_Fint* f_comm, const char* input_filename)
{
  // Initialize Kokkos if ELM hasn't already (for near future, it won't)
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
    ats_kokkos_init = true;
  }
  return reinterpret_cast<ELM_ATSDriver_ptr>(ATS::createELM_ATSDriver(f_comm, input_filename));
}


// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATSDriver_ptr ats)
{
  auto ats_ptr = reinterpret_cast<ATS::ELM_ATSDriver*>(ats);
  ats_ptr->finalize();
  delete ats_ptr;

  // If ATS initialized Kokkos, then finalize it here
  if (ats_kokkos_init && Kokkos::is_initialized()) {
    Kokkos::finalize();
    ats_kokkos_init = false;
  }
}


// call driver advance()
void ats_advance(ELM_ATSDriver_ptr ats,
                   double const * const dt,
                   bool const * const checkpoint,
                   bool const * const visualize)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(*dt, *checkpoint, *visualize);
}

// call driver advance_test()
void ats_advance_test(ELM_ATSDriver_ptr ats)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
}

void ats_get_mesh_info(ELM_ATSDriver_ptr ats,
                       int * const ncols_local,
                       int * const ncols_global,
                       int * const nlevgrnd)
{
  ATS::ELM_ATSDriver::MeshInfo info = reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->getMeshInfo();
  *ncols_local = info.ncols_local;
  *ncols_global = info.ncols_global;
  *nlevgrnd = info.nlevgrnd;
}


// call driver setup()
void ats_setup(ELM_ATSDriver_ptr ats)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setup();
}

// call driver initialize()
void ats_initialize(ELM_ATSDriver_ptr ats)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->initialize();
}

// call driver advance(dt)
void ats_advance(ELM_ATSDriver_ptr ats,
                   double dt,
                   bool checkpoint,
                   bool visualize)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(dt, checkpoint, visualize);
}

// call driver advance_test()
void ats_advance_test(ELM_ATSDriver_ptr ats)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(dt, checkpoint, visualize);
}


//
// Memory movement/data passing
// -----------------------------------------------------------------------------
double ats_get_scalar(ELM_ATSDriver_ptr ats, int scalar_id)
{
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->getScalar(reinterpret_cast<ATS::ELM::ScalarID>(scalar_id));
}

void ats_set_scalar(ELM_ATSDriver_ptr ats, int scalar_id, double in)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setScalar(reinterpret_cast<ATS::ELM::ScalarID>(scalar_id), in);
}

void ats_get_field(ELM_ATSDriver_ptr ats, int var_id, double * const in)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->getField(reinterpret_cast<ATS::ELM::VarID>(var_id), in);
}

double const * ats_get_field_ptr(ELM_ATSDriver_ptr ats, int var_id)
{
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->getFieldPtr(reinterpret_cast<ATS::ELM::VarID>(var_id));
}


double * ats_get_field_ptr_w(ELM_ATSDriver_ptr ats, int var_id)
{
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->getFieldPtrW(reinterpret_cast<ATS::ELM::VarID>(var_id));
}

void ats_set_field(ELM_ATSDriver_ptr ats, int var_id, double const * const in)
{
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->SetField(reinterpret_cast<ATS::ELM::VarID>(var_id), in);
}


#ifdef __cplusplus
}
#endif
