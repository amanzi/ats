
#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

#ifdef __cplusplus
extern "C"{
#endif

// allocate, call constructor and cast ptr to opaque ELM_ATS_DRIVER
ELM_ATS_DRIVER* ats_create() {
  return reinterpret_cast<ELM_ATS_DRIVER*>(new ATS::ELM_ATSDriver());
}
// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATS_DRIVER *ats) {
  delete reinterpret_cast<ATS::ELM_ATSDriver*>(ats);
}
// call driver setup()
void ats_setup(ELM_ATS_DRIVER *ats, MPI_Fint *f_comm, const char *input_filename) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setup(f_comm, input_filename);
}
// call driver initialize()
void ats_initialize(ELM_ATS_DRIVER *ats){
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->initialize();
}
// call driver advance(dt)
void ats_advance(ELM_ATS_DRIVER *ats, double *dt) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(dt);
}
// call driver advance_test()
void ats_advance_test(ELM_ATS_DRIVER *ats) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
}
#ifdef __cplusplus
}
#endif
