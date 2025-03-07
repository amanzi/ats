/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov), Xu Chonggang
*/

#ifndef PK_FATES_HH_
#define PK_FATES_HH_

#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "ISO_Fortran_binding.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include <string.h>

namespace Amanzi {
namespace BGC {

typedef struct {
  int nlevbed;
  int nlevdecomp;
  int patchno;
  int altmax_lastyear_indx_col;
  double temp_veg24_patch;
  double latdeg, londeg;
} site_info;

typedef struct {
  double dayl_factor; // scalar (0-1) for daylength
  double esat_tv;     // saturation vapor pressure at t_veg (Pa)
  double eair;        // vapor pressure of canopy air (Pa)
  double oair;        // Atmospheric O2 partial pressure (Pa)
  double cair;        // Atmospheric CO2 partial pressure (Pa)
  double rb;          // boundary layer resistance (s/m)
  double t_veg;       // vegetation temperature (Kelvin)
  double tgcm;        // air temperature at agcm reference height (Kelvin)
  double solad[2];    //direct radiation (W/m**2); 1=visible lights; 2=near infrared radition
  double solai[2];    //diffuse radiation (W/m**2); 1=visible lights; 2=near infrared radition
  double albgrd[2];   //!ground albedo (direct) 1=visiable; 2=near infrared (nir)
  double albgri[2];   //ground albedo (diffuse) 1=visiable; 2=near infrared (nir)
} PhotoSynthesisInput;


#ifdef __cplusplus
extern "C"
{
#endif
  void init_ats_fates(int*, site_info*);
  void init_soil_depths(int*, int*, site_info*, double*, double*, double*, double*);
  void init_coldstart(int*);
  void fatessetmasterproc(int*);
  void fatessetinputfiles(CFI_cdesc_t* clm, CFI_cdesc_t* fates);
  void fatesreadparameters();
  void fatesreadpfts();
  void set_fates_global_elements();
  void get_nlevsclass(int*);
  void dynamics_driv_per_site(int*,
                              int*,
                              site_info*,
                              double*,
                              double*,
                              double*,
                              double*,
                              double*,
                              double*);

  void wrap_btran(int*, double*, double*, double*, double*, double*);
  void wrap_photosynthesis(double*, double*, int*, double*, PhotoSynthesisInput*);
  void wrap_sunfrac(int* array_size, double* forc_solad, double* forc_solai);
  void wrap_canopy_radiation(double* jday, int* array_size, double* albgrd, double* albgri);

  void calculate_biomass(double* ats_biomass_array, int nsites, int num_scls);

#ifdef __cplusplus
} // extern "C"
#endif

class FATES_PK : public PK_Physical {
 public:
  FATES_PK(Teuchos::ParameterList& FElist,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~FATES_PK() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  //virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  //virtual void initialize(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);


  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S){};

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // -- provide a timestep size
  virtual double get_dt();

  virtual void set_dt(double dt) { dt_ = dt; }

 protected:
  void
  FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec, double* col_vec, int ncol);
  void ColDepthDz_(AmanziMesh::Entity_ID col,
                   Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                   Teuchos::Ptr<Epetra_SerialDenseVector> dz);


  double dt_, dt_photosynthesis_, dt_site_dym_;
  double t_photosynthesis_, t_site_dym_;

  bool surface_only_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_, mesh_domain_;
  Key domain_surf_;
  Key trans_key_;
  Key precip_key_, air_temp_key_, vp_air_key_, wind_key_, co2a_key_;
  Key poro_key_, sat_key_, suc_key_, soil_temp_key_;
  Key met_decomp_key_, cel_decomp_key_, lig_decomp_key_;
  Key longwave_key_, incident_rad_key_;

  std::vector<double> t_soil_;   // soil temperature
  std::vector<double> vsm_;      // volumetric soil moisture vsm_ = S * poro;
  std::vector<double> poro_;     // porosity
  std::vector<double> eff_poro_; //effective porosity  = porosity - vol_ice
  std::vector<double> suc_;      //suction head

  int patchno_, nlevdecomp_, nlevsclass_;
  int ncells_owned_, ncells_per_col_, clump_;
  std::vector<site_info> site_;

  // factory registration
  static RegisteredPKFactory<FATES_PK> reg_;
};

} // namespace BGC
} // namespace Amanzi

#endif
