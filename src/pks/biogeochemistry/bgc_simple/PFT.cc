/*
  Plant Functionality Type class

  -- copy constructor of array?
  -- Epetra_Vector copy constructor copy values?

*/

#include "utils.hh"

#include "PFT.hh"

namespace Amanzi {
namespace BGC {

PFT::PFT(std::string pft_type_, int ncells) :
    pft_type(pft_type_),
    BRootSoil(ncells) {}

PFT::PFT(std::string pft_type_, int ncells, double* brootcells_) :
    pft_type(pft_type_),
    BRootSoil(View, brootcells_, ncells) {}

void PFT::Init(double col_area)
{
  Bleaf = 0.;
  Bleafmemory = 0.;
  Broot = 0.;
  Bstem = 0.;
  Bstore = 0.;

  
  SLA = 30;
  leaflitterfrc[0] = 0.1;
  leaflitterfrc[1] = 0.5;
  leaflitterfrc[1] = 0.4;
  rootlitterfrc[0] = 0.1;
  rootlitterfrc[1] = 0.5;
  rootlitterfrc[1] = 0.4;
  stemlitterfrc[0] = 0.05;
  stemlitterfrc[1] = 0.15;
  stemlitterfrc[1] = 0.8;

  leaf2rootratio = 1.0;
  leaf2stemratio = 5;
  root2leafrespratio = 0.8;
  stem2leafrespratio = 0.05;
  rootlongevity = 10.0;
  leaflongevity = 4.0;
  stemlongevity = 100.0;

  Vcmax25 = 100;
  Emax25 = 0.1*Vcmax25;
  Jmax25 = 200;
  
  GDDleafon = 100;
  GDDbase = 0.0;
  GDD = 0.0;

  mResp = 0.;
  gResp = 0.;

  annNPP = 0.0;
  GPP = 0.0;
  NPP = 0.0;
  ET = 0.0; 

  leafstatus = 1;
  evergreen = false;
  
  LUE = 0.06;
  LER = 0.8;
  mp = 9.0;
  swpo = -0.5;
  swpc = -2.5; 

  lai = 0.;
  laimemory = 0.;
  totalBiomass = 0.;
  maxRootD = 0.5;
  rootD = 0.;

  minLeafWP = -2.0;
  storagecleaf2sw = 1.5;
  storagecroot2sw = 15;
  tar_leafstorageccon = 0.2;
  gRespF = 0.25;
  storagecRspFrc = 0.8;

  bleafon = 0.;
  bleafoff = 0;
  leafondays = 5;
  leafoffdays = 5;
  leafoffdaysi = 365.25;
  leafondaysi = 0;

  CSinkLimit = 0.0;
  seedrainlai = 0.01;
  maxLAI = 1;

  for (int i = 0; i < 10; i++){
    annCBalance[i] = 0;
  }
  for (int c = 0; c != BRootSoil.Length(); ++c) {
    BRootSoil[c] = 0.;
  }
}

void PFT::Init(Teuchos::ParameterList& plist,double col_area)
{
  Init(col_area);
  // ex
  // rootlongevity = plist.get<double>("root longevity", 4.0);
  //note default vals below are those of sedge
  maxRootD = plist.get<double>("max root depth [m]", 0.5);
  Vcmax25 = plist.get<double>("Vcmax25 [micromol CO2 / m^2(leaf) s]", 100.);
  Jmax25 = 2.0 * Vcmax25; 
  Emax25 = plist.get<double>("Emax25 [micromol C / m^2(leaf) s]", 10.);
  SLA = plist.get<double>("SLA [m^2(leaf) / kg C", 16);
  evergreen = plist.get<bool>("evergreen", false);
  leaf2rootratio =  plist.get<double>("ratio of leaf to root [-]", 1.0);
  leaf2stemratio =  plist.get<double>("ratio of leaf to stem [-]", 5.0);

  leaflongevity =  plist.get<double>("leaf longevity [yr]", 4.0);
  rootlongevity =  plist.get<double>("root longevity [yr]", 4.0);
  stemlongevity =  plist.get<double>("stem longevity [yr]", 10.0);

  mp =  plist.get<double>("stomatal conductance to photosynthetic rate slope (mp) [?]", 9.0);
  

  Bleaf = 1.0/SLA * col_area;//es note that all the following B vals are per m^2, so that elsewhere they should be *gridarea to account for varying grid areas.
  Bstore = 2* Bleaf;
  Bleafmemory = Bleaf;
  Bstem = Bleaf/leaf2stemratio;
  Broot = Bleaf/leaf2rootratio;
  totalBiomass = Bleaf + Bstem + Broot;

  if (evergreen) { 
    Bleafmemory = 0.;
  }else{
    Bleaf = 0.;
  }

}


// Initialize the root distribution
void PFT::InitRoots(const Epetra_SerialDenseVector& SoilTArr,
                    const Epetra_SerialDenseVector& SoilDArr,
                    const Epetra_SerialDenseVector& SoilThicknessArr) {

  // locate the root depth
  int nSoilLayers = SoilTArr.Length();
  double initPermD = PermafrostDepth(SoilTArr,SoilThicknessArr,273.15);
  //rootD = std::min(maxRootD, initPermD);
  rootD = 0.8*maxRootD; //es following Xu's advice - 0.8 * maxRootD accounts for the previous year root depth that will contribute to this yr Broot. Note that if spinning data exists feed in last year's root depth instead of this.
 
  double totalweights = 0.0;
  for (int c=0; c!=nSoilLayers; ++c) {
    if (SoilDArr[c] < rootD){
      BRootSoil[c] = SoilThicknessArr[c]*std::exp((rootD - SoilDArr[c]) / rootD);
      totalweights += BRootSoil[c];
    } else {
      BRootSoil[c] = 0.0;
    }
  }

  if (BRootSoil[0] <= 0.0) { //avoid the numerical errors
    BRootSoil[0] = 1.0;
    totalweights = 1.0;
  }

  for (int c=0; c!=nSoilLayers; ++c) {
    BRootSoil[c] = BRootSoil[c] / totalweights * Broot;
  }

  AssertRootBalance_or_die();
  return;
}

} // namespace
} // namespace
