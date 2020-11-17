/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/
//! Plant Functionality Type class
/*!

Paramters for a plant functional type.

.. _pft-spec:
.. admonition:: pft-spec

    * `"max root depth [m]`" ``[double]`` Depth of the roots.
    * `"Vcmax25 [micromol CO2 / m^2(leaf) s]`" ``[double]``
    * `"Emax25 [micromol C / m^2(leaf) s]`" ``[double]``
    * `"SLA [m^2(leaf) / kg C`" ``[double]`` Specific leaf area.
    * `"evergreen`" ``[bool]`` **false** Is evergreen?
    * `"ratio of leaf to root [-]`" ``[double]``
    * `"ratio of leaf to stem [-]`" ``[double]``
    * `"leaf longevity [yr]`" ``[double]``
    * `"root longevity [yr]`" ``[double]``
    * `"stem longevity [yr]`" ``[double]``

    * `"height [m]`" ``[double]``
    * `"stomatal conductance to photosynthetic rate slope (mp) [?]`" ``[double]``
    * `"wood density [kg m^-3]`" ``[double]``

*/

#ifndef ATS_BGC_PFT_HH_
#define ATS_BGC_PFT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"

#include "dbc.hh"

namespace Amanzi {
namespace BGC {

class PFT {

 public:
  PFT(std::string pft_type, int ncells);
  PFT(std::string pft_type, int ncells, double* brootsoil_);

  virtual Teuchos::RCP<PFT> Clone() const {
    return Teuchos::rcp(new PFT(*this));
  }
  virtual ~PFT() {}

  // Default copy constructor does the right thing.
  // Default operator= does the right thing.

  virtual void Init(double col_area);
  void Init(Teuchos::ParameterList& plist, double col_area);

  void InitRoots(const Epetra_SerialDenseVector& SoilTArr,
                 const Epetra_SerialDenseVector& SoilDArr,
                 const Epetra_SerialDenseVector& SoilThicknessArr);

  bool AssertRootBalance_or_die() {
    double totalRootW = BRootSoil.Norm1();
    AMANZI_ASSERT(std::abs(totalRootW - Broot) < 1.e-6);
    return std::abs(totalRootW - Broot) < 1.e-6;
  }


  double Bleaf; //kg c/m2 land
  double Bleafmemory; //kg c/m2 land, leaf biomass for the last year
  double Broot; //kg c/m2 land
  double Bstem; //kg c/m2 land
  double Bstore; //kg c/m2 land, the storage carbon, not functional yet.
  double SLA; //leaf mass per unit area (m2 leaf/kg c)
  double leaflitterfrc[3]; //[0],[1],[2] leaf suguar, celluose, and lignin fractions
  double stemlitterfrc[3]; //[0],[1],[2] stem suguar, celluose, and lignin fractions
  double rootlitterfrc[3]; //[0],[1],[2] root suguar, celluose, and lignin fractions
  double leaf2rootratio; //leaf to root biomass ratio
  double leaf2stemratio; //leaf to stem biomass ratio
  double root2leafrespratio; //leaf to root respiration ratio
  double stem2leafrespratio; //leaf to stem respiration ratio
  double rootlongevity; //years
  double leaflongevity; //years
  double stemlongevity; //years
  double Vcmax25; //umol CO2/m2 leaf/s
  double Emax25; ///export rate of carbon (umol c export/s/m2 leaf) at 25 degrees
  double Jmax25; //umol electron/m2 leaf/s
  double GDDleafon; //growing degree days for leaf on
  double GDD; //growing degree days
  double mResp; //maintanence respiration, kg C/grid cell /time step
  double gResp; //maintanence respiration, kg C/grid cell /time step
  double GDDbase; //growing degree days base temperature
  double annNPP;
  double GPP;//gross primary proeduction; kg C/grid cell /time step
  double NPP;//net primary proeduction; kg C/grid cell/time step
  double ET; //Transpiration; kg water/time step/grid cell
  int leafstatus; //1=off; 2==on
  bool evergreen;  //1=yes; 0=no
  double LUE; //light use efficiency
  double LER; //light exptinction rate
  double mp; //the slope between stamta conductance and photosynthetic rate
  double swpo; //the soil water potential plant fully open stomata (MPa)
  double swpc; //the soil water potential plant fully close stomata (MPa)
  double lai; //leaf are index
  double laimemory; //memorized leaf area index
  double totalBiomass; //alive biomass;
  double maxRootD; //maximum root depth (m)
  double rootD; //curent root depth(m)
  double minLeafWP; //minimum leaf water potential
  double storagecleaf2sw; //the ratio of storage carbon  in leaves compared to stem, given the same amount of drymass of leaf and stem (exluding storage)
  double storagecroot2sw; //the ratio of storage carbon in  root compared to stem, given the same amount of drymass of root and stem (exluding storage)
  double tar_leafstorageccon; // targe storage carbon concentration in leaves, g storage c / g leaf c (including storage)
  double gRespF; //growth respiration fraction
  double storagecRspFrc; // the fraction of storage carbon that are recycled
  double bleafon;
  int leafondays; //amount of leaf on per day
  int leafoffdays;
  int bleafoff; //amount of leaf fall per day
  double leafondaysi; //count of the number of days leaves on
  double leafoffdaysi; //cound the number of days leaves off
  double CSinkLimit; //the carbon sink limitation
  double seedrainlai; //the seed rain in cases of dying off
  int maxLAI; //maximum leaf area index;
  double annCBalance[10];//annual carbon balance for each leave layers
  double height;
  double carbon2biomass; // 2, the factor to convert Kg carbon to Kg biomass
  double wood_density;
  double aroot_radius; // m fine root radius
  
  Epetra_SerialDenseVector BRootSoil; //root distribution at different soil depths
  std::string pft_type;

};


} // namespace
} // namespace

#endif

