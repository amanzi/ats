
#ifdef __cplusplus
extern "C" {
    class ELM_ATSDriver;
    typedef ELM_ATSDriver ELM_ATS_DRIVER;
#else
    typedef struct ELM_ATS_DRIVER ELM_ATS_DRIVER;
#endif

// allocate, call constructor and cast ptr to opaque ELM_ATS_DRIVER
ELM_ATS_DRIVER* ats_create();
// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATS_DRIVER *ats);
// call driver setup()
void ats_setup(ELM_ATS_DRIVER *ats, MPI_Fint *f_comm, const char *input_filename);
// call driver initialize()
void ats_initialize(ELM_ATS_DRIVER *ats);
// call driver advance(dt)
void ats_advance(ELM_ATS_DRIVER *ats, double *dt);
// call driver advance_test()
void ats_advance_test(ELM_ATS_DRIVER *ats);
// call driver set_sources()
void ats_set_sources(ELM_ATS_DRIVER *ats, double *soil_infiltration, double *soil_evaporation,
    double *root_transpiration, int *ncols, int *ncells);

#ifdef __cplusplus
}
#endif
