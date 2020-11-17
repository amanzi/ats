/*
  Main functions for biogeochemistry on a column.

  Author: Ethan Coon (ecoon@lanl.gov)
  Chonggang Xu (cxu@lanl.gov)


  Issues:

  -- has no knowledge of snow cover?
  -- has no use of wind speed ref ht

  License: BSD
*/

#include <algorithm>
#include <cmath>

#include "vegetation.hh"
#include "bgc_simple_funcs.hh"

namespace Amanzi {
namespace BGC {

  // t and dt [s]
  // gridarea [m^2]
  // qSWin [W/m^2]
  // tair [K]
  // met.windv [m/s]
  // met.relhum [-]
  // met.CO2a [ppm]
  // SoilTArr [K] (soil temperature)
  // SoilWPArr [Pa] (water pressure)
  // SoilDArr [m] (depth)
  // SoilThicknessArr [m] (dz)
  // TransArr[kg H2O/m3/s]
  // sw_shaded[W/m^s] (shaded shortwave radiation that makes it to the surface)
  void BGCAdvance(double t, double dt, double gridarea, double cryoturbation_coef,
		  const MetData& met,
		  const Epetra_SerialDenseVector& SoilTArr,
		  const Epetra_SerialDenseVector& SoilWPArr,
		  const Epetra_SerialDenseVector& SoilDArr,
		  const Epetra_SerialDenseVector& SoilThicknessArr,
		  std::vector<Teuchos::RCP<PFT> >& pftarr,
		  std::vector<Teuchos::RCP<SoilCarbon> >& soilcarr,
		  Epetra_SerialDenseVector& SoilCO2Arr,
		  Epetra_SerialDenseVector& TransArr,
		  double& sw_shaded)
  {
  // required constants
  double p_atm = 101325.;
  double Cv = 1.2e-8;
  double dt_days = dt / 86400.;
  double t_days = t / 86400.;
  double wp_max = -1.e-6; //MPa wp = p - p_atm
  double wp_min = -10.; //MPa
  int max_leaf_layers = 10;
  int ncells = SoilTArr.Length();

  // calculate fractional day length
  int doy = std::floor(std::fmod(t_days, 365.25));
  if (doy == 0) doy = 365;
  double daylen = DayLength(met.lat, doy);

  // initialize transpiration to zero
  for (int k=0; k!=ncells; ++k) TransArr[k] = 0.0 ;

  // convert to PAR at the daytime
  const double PAR = met.qSWin * 2.3 * 24.0 *60.0 / daylen;

  // determine the thaw depth
  const double thawD = PermafrostDepth(SoilTArr,SoilThicknessArr,273.15);


  //----------------------------------------------------------------------
  // loop through the list of PFTs

  for (std::vector<Teuchos::RCP<PFT> >::iterator pft_iter=pftarr.begin();
       pft_iter!=pftarr.end(); ++pft_iter) {

    PFT& pft = *(*pft_iter);
    pft.GPP = 0.0;
    pft.ET = 0.0;

    if (pft.totalBiomass > 0.0) {
      if (!pft.evergreen) {
        //---------------------------------------------------------------------------------
        //calculate plant phenology
        if (pft.leafstatus == 1) {
          if (met.tair - 273.15 > pft.GDDbase) {
            pft.GDD += (met.tair - 273.15 - pft.GDDbase) * dt_days;
          }

          if (pft.GDD >= pft.GDDleafon){
            pft.leafstatus = 2;
            pft.bleafon = std::min(pft.Bleafmemory,
                    pft.Bstore - 0.015*pft.Bleafmemory) / pft.leafondays;

            if (pft.bleafon < 0.0)
              pft.bleafon = 0.5* pft.Bstore / pft.leafondays;

            pft.leafondaysi = 0;
            pft.Bleafmemory = pft.bleafon*pft.leafondays;
          }

          if (pft.leafoffdaysi < pft.leafoffdays) {
            pft.leafoffdaysi += dt_days;
          }
        }

        if (pft.leafstatus == 2) {
          // add leaves
          if (pft.leafondaysi < pft.leafondays){
            double bleafon = std::min(pft.bleafon, 0.9*pft.Bstore);
            pft.Bleaf += bleafon * dt_days;
            pft.leafondaysi += dt_days;
            pft.Bleafmemory -= bleafon * dt_days;
            pft.Bstore -= bleafon * dt_days;
          }

          if (pft.leafondaysi >= pft.leafondays && pft.leafondaysi < 365.25) {
            pft.Bleafmemory = 0.0;
            pft.bleafon = 0.0;

            // AFFECTED BY DT! --etc
            pft.leafondaysi = 365.25;
          }

          // shed leaves
          // AFFECTED BY DT! --etc
          if (daylen < 655 && pft.leafoffdaysi >= 365.25) {
            pft.bleafoff = pft.Bleaf / pft.leafoffdays;
            pft.leafoffdaysi = 0;
          }

          // shed leaves, gradual leaf fall
          if (pft.leafoffdaysi < pft.leafoffdays) {
            pft.Bleaf -= pft.bleafoff * dt_days;
            pft.Bleafmemory += pft.bleafoff * dt_days;

            // litter transfer
            double stemstoragecratio = pft.Bstore
                / ((pft.Bleaf + pft.Bleafmemory)*pft.storagecleaf2sw + pft.Bstem + pft.Broot*pft.storagecroot2sw);
            double leafstoragecratio = pft.storagecleaf2sw*stemstoragecratio;
            double storagecdrawn = pft.bleafoff * dt_days * leafstoragecratio*(1.0 - pft.storagecRspFrc);
            pft.Bstore -= storagecdrawn;
            soilcarr[0]->SOM[0] += storagecdrawn;
            double carbondrawnLeaf = pft.bleafoff * dt_days;

            for (int k=0; k!=2; ++k) {
              soilcarr[0]->SOM[k+1] += carbondrawnLeaf*pft.leaflitterfrc[k];
            }

            pft.leafoffdaysi += dt_days;
          }

          // shed leaves complete, final check
          if (pft.leafoffdaysi >= pft.leafoffdays && pft.leafoffdaysi < 365.25) {
            pft.leafstatus = 1;
            pft.Bleafmemory += pft.Bleaf * dt_days;
            pft.leafoffdaysi = 365.25;
            pft.GDD = 0.0;

            // litter transfer
            double stemstoragecratio = pft.Bstore
                / ((pft.Bleaf + pft.Bleafmemory)*pft.storagecleaf2sw + pft.Bstem + pft.Broot*pft.storagecroot2sw);
            double leafstoragecratio = pft.storagecleaf2sw*stemstoragecratio;
            double storagecdrawn = pft.Bleaf * leafstoragecratio * (1.0 - pft.storagecRspFrc);
            pft.Bstore -= storagecdrawn;
            soilcarr[0]->SOM[0] += storagecdrawn;
            double carbondrawnLeaf = pft.Bleaf;

            for (int k=0; k!=2; ++k) {
              soilcarr[0]->SOM[k+1] += carbondrawnLeaf*pft.leaflitterfrc[k];
            }
            pft.Bleaf = 0.0;
          }
        }
      }

      pft.lai = pft.Bleaf*pft.SLA / gridarea;
      pft.laimemory = pft.Bleafmemory*pft.SLA / gridarea;
      pft.totalBiomass = pft.Bleaf + pft.Broot + pft.Bstem + pft.Bstore;

      //---------------------------------------------------------------------------------
      //calculate leaf projections and light attenuations
      double PARi0 = PAR;
      // note loop dimensions range from 0 to current PFT.  This could be
      // tracked and refactored out. --etc
      for (std::vector<Teuchos::RCP<PFT> >::iterator pft_other_iter=pftarr.begin();
           pft_other_iter!=pft_iter; ++pft_other_iter) {
        PARi0 *= std::exp(-(*pft_other_iter)->LER * (*pft_other_iter)->lai);
      }
      double PARi = PARi0;
      //----------------------------------------------------------------
      // photosynthesis and respiration, no soil water limitation yet and need
      // to be implemented for more realistic simulation
      const double nleaflayers = std::ceil(pft.lai);
      double leafresptotal = 0.0;

      double psn = 0.;
      double tleaf = 0.;
      double leafresp = 0.;
      double ET = 0.; //umol H2O/m2 leaf/s
      pft.GPP = 0.;
      pft.ET = 0.;
      //=============================================================================
      //distribution water limitations
      double Btran = 0.0;
      double soil_wp;
      double WFactor;
      double rootFrac;
      for (int k=0; k!=ncells; ++k) {
            soil_wp = std::max(std::min( (SoilWPArr[k] - p_atm) / 1.e6, pft.swpo), pft.swpc);
	    WFactor = (pft.swpc-soil_wp)/(pft.swpc-pft.swpo);
	    rootFrac = pft.BRootSoil[k] / pft.Broot;
            Btran = Btran + WFactor * rootFrac;
      }
      //=========================================================================
      //do photosynthesis and respirations
      double  Vcmax25 = pft.Vcmax25 * (1.0 - pft.CSinkLimit);
      double  relRad; // relative radiation compared to the top of canopy
      double  relCLNCa; // relative leaf nitrogen content comapred to the top canopy
      double  Vcmax25i; // Vcmax at each leaf layers 
      if (PAR>0.0) {
	  for (int leaf_layer=0; leaf_layer!=max_leaf_layers; ++leaf_layer){
	    relRad   = PARi/PAR;
	    relCLNCa = 0.1802 * std::log(relRad)+1.0; //see Ali et al 2015
	    relCLNCa = std::max(0.2,relCLNCa);
	    relCLNCa = std:: min(1.0,relCLNCa);
	    Vcmax25i = Vcmax25 * relCLNCa;
	    Photosynthesis(PARi, pft.LUE, pft.LER, p_atm,
			   met.windv,double(met.tair - 273.15), met.relhum, met.CO2a, pft.mp, Vcmax25i,
			   &psn, &tleaf, &leafresp,&ET);
	    psn *= Btran;
	    ET  *= Btran; 
	    if (thawD <= 0.0) {
		psn = 0.0;
		ET = 0.0;
	    }         
	    leafresp *= dt * Cv;
	    double rootresp = leafresp / pft.leaf2rootratio * pft.root2leafrespratio;
	    double stemresp = leafresp / pft.leaf2stemratio * pft.stem2leafrespratio;
	    double NPP = daylen * 60.0 * Cv * psn * dt_days - leafresp - rootresp - stemresp;
	    pft.annCBalance[leaf_layer] += NPP;
	    if (leaf_layer < nleaflayers) {
	      if (leaf_layer == nleaflayers - 1){
		psn *= (pft.lai - nleaflayers + 1);
		leafresp *= (pft.lai - nleaflayers + 1);
		ET *= (pft.lai - nleaflayers + 1);
	      }
	      pft.GPP = pft.GPP +  daylen * 60.0 * Cv* psn * gridarea * dt_days;
	      leafresptotal += dt * Cv * leafresp * gridarea;
	      pft.ET += dt * 18.0 * 1.0e-9 * ET * gridarea;
	    }
	    PARi *= std::exp(-pft.LER);
         }
      }
      //------------------------------------------------------------------------------------------------
      // respiration
      double stemresp = 0.;
      double rootresp = 0.;

      if (pft.GPP > 0.0){
        stemresp = leafresptotal*pft.Bstem / pft.Bleaf*pft.stem2leafrespratio;
        rootresp = 0.0;
        double refTFactor = TEffectsQ10(2.0, tleaf, 25.0);
        for (int k=0; k!=ncells; ++k){
          double TFactor = TEffectsQ10(2.0, SoilTArr[k] - 273.15, 25.0);
          rootresp += leafresptotal*pft.BRootSoil[k] / pft.Bleaf * TFactor / refTFactor
              * pft.root2leafrespratio;
        }

      } else {
        PARi = PARi0;
        Vcmax25i = Vcmax25;

        Photosynthesis(PARi, pft.LUE, pft.LER, p_atm,
                       met.windv,double(met.tair - 273.15), met.relhum, met.CO2a, pft.mp, Vcmax25i,
                       &psn, &tleaf, &leafresp, &ET);

        if(met.tair<273.15) leafresp = leafresp/10.0; //winter hypbernation

        leafresp *= dt * Cv;
        rootresp = leafresp / pft.leaf2rootratio * pft.root2leafrespratio;
        stemresp = leafresp / pft.leaf2stemratio * pft.stem2leafrespratio;
        double NPP = - leafresp - rootresp - stemresp;

        for (int leaf_layer=0; leaf_layer!=max_leaf_layers; ++leaf_layer) {
          relCLNCa = -0.1802 * pft.LER * leaf_layer +1.0; //see Ali et al 2015
          relCLNCa = std::max(0.2,relCLNCa);
          relCLNCa = std:: min(1.0,relCLNCa);
          pft.annCBalance[leaf_layer] += NPP *  relCLNCa ;

        }

        leafresptotal = dt * Cv*leafresp*gridarea;

        double bleaf0 = 1.0 / pft.SLA * gridarea;
        stemresp = leafresptotal*pft.Bstem / bleaf0*pft.stem2leafrespratio;
        rootresp = 0.0;
        double refTFactor = TEffectsQ10(2.0, tleaf, 25.0);
        for (int k=0; k!=ncells; ++k) {
          double TFactor = TEffectsQ10(2.0, SoilTArr[k] - 273.15, 25.0);
          rootresp += leafresptotal*pft.BRootSoil[k] / bleaf0 * TFactor / refTFactor
              * pft.root2leafrespratio;
        }
        leafresptotal = 0.0;
      }

      pft.mResp = leafresptotal + stemresp + rootresp;

      //------------------------------------------------------------------------
      // calculate plant allocations, biomass allocations to root, leaf, stem and storage
      if (pft.mResp > 0.1*pft.Bstore) {
        //avoid negative carbon fluxes and downregulate the maintenance respiration
        pft.mResp = 0.05*pft.Bstore;
      }

      pft.Bstore = pft.Bstore - pft.mResp;
      double stemstoragecratio = pft.Bstore / ((pft.Bleaf + pft.Bleafmemory)
              * pft.storagecleaf2sw + pft.Bstem + pft.Broot*pft.storagecroot2sw);
      double leafstoragecratio = pft.storagecleaf2sw*stemstoragecratio;

      //conversion from ratio to concentration
      double leafstorageccon = leafstoragecratio > 0. ? leafstoragecratio / (1.0 + leafstoragecratio) : 0.;

      //-------------------------------------------------------------------------------------- -
      // carbon sink rate based on carbon storage
      pft.Bstore = pft.Bstore + pft.GPP;
      double frac = leafstorageccon / pft.tar_leafstorageccon;
      double csink_factor = 1.0 - frac > 0. ? std::exp(-std::pow(frac, 3.0)) : 1.0;

      double Emax = met.tair - 273.15 > 0. ? TEffectsQ10(2.0, met.tair - 273.15, 25.0) * pft.Emax25 : 0.;

      // calculate the carbon sink limitation for photosynthesis
      pft.CSinkLimit = 0.0;
      if (leafstorageccon>2 * pft.tar_leafstorageccon){
        pft.CSinkLimit = (leafstorageccon - 2.0 * pft.tar_leafstorageccon)
            / pft.tar_leafstorageccon;
        pft.CSinkLimit = std::min(1.0, pft.CSinkLimit);
      }

      double GrowthFlux = 0.0010368* Emax*csink_factor*pft.lai*gridarea*dt_days;
      // -- 0.0010368 is a conversion factor from umol C/m2/s->kg C/m2/day
      // for carbon storage boundary check,avoid numerical errors
      GrowthFlux = std::min(GrowthFlux, 0.5*pft.Bstore);

      pft.Bstore = pft.Bstore - GrowthFlux;
      pft.gResp = pft.gRespF* GrowthFlux;
      GrowthFlux = (1.0 - pft.gRespF)*GrowthFlux;

      //-------------------------------------------------------------------------------------------------------
      // carbon allocations
      if (GrowthFlux > 0.0) {
        double grwBleaf = 0.0;
        double grwBroot = 0.0;
        double grwBstem = 0.0;
        double totalNonStoreB = pft.Bleaf + pft.Bleafmemory + pft.Broot + pft.Bstem;
        double tarBleaf = totalNonStoreB
            / (1.0 / pft.leaf2rootratio + 1.0 + 1.0 / pft.leaf2stemratio);

        double tarLAI = tarBleaf*pft.SLA/gridarea;
        if (tarLAI > pft.maxLAI){
          tarBleaf = pft.maxLAI / pft.SLA*gridarea;
        }
        double tarBroot = tarBleaf / pft.leaf2rootratio;
        double tarBstem = tarBleaf / pft.leaf2stemratio;
        double tarBtotal = tarBleaf + tarBstem + tarBroot;

        double deficitBleaf = std::max(tarBleaf - pft.Bleaf - pft.Bleafmemory, 0.);
        double deficitBroot = std::max(tarBroot - pft.Broot, 0.);
        double deficitBstem = std::max(tarBstem - pft.Bstem, 0.);
        double totalDeficit = deficitBstem + deficitBroot + deficitBleaf;

        double frcBleaf = totalDeficit > 0 ? deficitBleaf / totalDeficit : 0.;
        double frcBroot = totalDeficit > 0 ? deficitBroot / totalDeficit : 0.;
        double frcBstem = totalDeficit > 0 ? deficitBstem / totalDeficit : 0.;

        if (totalDeficit >= GrowthFlux){
          grwBleaf += frcBleaf*GrowthFlux;
          grwBroot += frcBroot*GrowthFlux;
          grwBstem += frcBstem*GrowthFlux;
          GrowthFlux = 0.0;
        } else {
          grwBleaf += deficitBleaf;
          grwBroot += deficitBroot;
          grwBstem += deficitBstem;
          GrowthFlux -= totalDeficit;

          //--------------------------------------------------------------------
          //grow
          if (pft.lai <= pft.maxLAI){
            frcBleaf = tarBtotal > 0 ? tarBleaf / tarBtotal : 0.;
            frcBroot = tarBtotal > 0 ? tarBroot / tarBtotal : 0.;
            frcBstem = tarBtotal > 0 ? tarBstem / tarBtotal : 0.;

            grwBleaf += frcBleaf*GrowthFlux;
            grwBroot += frcBroot*GrowthFlux;
            grwBstem += frcBstem*GrowthFlux;
          } else {
            // put back to storage
            pft.Bstore = pft.Bstore + GrowthFlux / (1 - pft.gRespF);
            pft.gResp = pft.gResp - pft.gRespF*GrowthFlux / (1 - pft.gRespF);
          }

          pft.Bleaf = pft.Bleaf + grwBleaf;
          pft.Bstem = pft.Bstem + grwBstem;
        }

        pft.NPP = pft.GPP - pft.gResp - pft.mResp;
        pft.annNPP = pft.annNPP + pft.NPP;

        //================================================================================
        //calculate new root distribution based on growth
        //caculate the root depth (90% root)
        if (grwBroot > 0.0){
          double totalweights = 0.0;
          bool findflag = false;

          for (int k=0; k!=ncells && !findflag; ++k) {
            totalweights += pft.BRootSoil[k] / pft.Broot;
            if (totalweights >= 0.85) {
              pft.rootD = SoilDArr[k];
              findflag = true;
            }
          }
          AMANZI_ASSERT(findflag);

          //-------------------------------------------------------------
          // calculate the root growth distribution in soil layers, the
          // distribution is based growth rate, water availability and
          // currrent root biomass

          totalweights = 0.0;
          std::vector<double> weightArr(ncells, 0.);

          for (int k=0; k!=ncells; ++k) {
            double TFactor = SoilTArr[k] < 273.15 ? 0. :
                TEffectsQ10(2.0, SoilTArr[k] - 273.15, 25.0) * HighTLim(SoilTArr[k]-273.15);
            double soil_wp = std::max(std::min( (SoilWPArr[k] - p_atm) / 1.e6, wp_max), wp_min);
            double WFactor = std::max((soil_wp - pft.minLeafWP) / (-0.05 - pft.minLeafWP),0.);
            if (SoilDArr[k] >= thawD) {
	      WFactor = 0.0;
            }
            weightArr[k] = pft.BRootSoil[k] * TFactor*WFactor;
            totalweights += weightArr[k];
          }

          if (totalweights <= 0.0){ //completely frozen
            totalweights = pft.Broot;
            for (int k=0; k!=ncells; ++k) {
              weightArr[k] = pft.BRootSoil[k];
            }
          }

          for (int k=0; k!=ncells; ++k) {
            weightArr[k] /= totalweights;
          }
          //---------------------------------------
          // Check root mass balance
          double sumWeight = 0.0 ;
          for (int k=0; k!=(ncells - 1); ++k) {
            sumWeight = sumWeight + weightArr[k];
          }
          pft.AssertRootBalance_or_die();

          //-------------------------------------------------------------
          // determine the growth direction depending on the rooting depth
          // grow downwards
          if (pft.rootD < thawD && pft.rootD < pft.maxRootD) {
            for (int k=0; k!=(ncells - 1); ++k) {
              if (SoilDArr[k] < thawD) {
                pft.BRootSoil[k+1] = pft.BRootSoil[k+1] + grwBroot*weightArr[k];
              }
            }
            if (SoilDArr[ncells-1] < thawD) {
              //bottom soil layer, stay put
              pft.BRootSoil[ncells-1] += grwBroot*weightArr[ncells-1];
            }
          } else {
            // grow horizontally
            for (int k=0; k!=ncells; ++k) {
              pft.BRootSoil[k] += grwBroot*weightArr[k];
            }
          }

          pft.Broot = pft.Broot + grwBroot;

          // Check root mass balance
          pft.AssertRootBalance_or_die();
        }
      }

      //=====================================================================================================
      // vegetation mortality and tissue turn over
      double mort = 0.0;
      // do plant mortality---very simple approach, but may improve later for more mechanistic approaches
      if (pft.Bstore < 0.01*(pft.Bleaf + pft.Bleafmemory)) {
        // kill vegetation
        mort = dt_days * 0.1 / 365.25; // AFFECTED BY DT! --etc
      }

      if ((pft.lai + pft.laimemory) < 0.001 ||
          pft.Bstore < 0.00001*(pft.Bleaf + pft.Bleafmemory)) {
        // kill all to avoid very small vegetation types and numerical errors
        mort = 1.0;
        std::cout << "WARNING: plant killed for pft " << pft.pft_type;
      }

      if ( mort > 0.0) {
        //-----------------------------------------------------------------------------------------------------
        // transfer to litter, due to mortality
        double stemstoragecratio = pft.Bstore / ((pft.Bleaf + pft.Bleafmemory)*pft.storagecleaf2sw + pft.Bstem + pft.Broot*pft.storagecroot2sw);
        double leafstoragecratio = pft.storagecleaf2sw*stemstoragecratio;
        double rootstoragecratio = pft.storagecroot2sw*stemstoragecratio;
        double carbondrawnLeaf = mort*pft.Bleaf;
        double carbondrawnStem = mort*pft.Bstem;

        soilcarr[0]->SOM[0] += carbondrawnLeaf*leafstoragecratio;
        //pft.Bstore = pft.Bstore - carbondrawnLeaf*leafstoragecratio;

        soilcarr[0]->SOM[0] += carbondrawnStem*stemstoragecratio;
        //pft.Bstore = pft.Bstore - carbondrawnStem*stemstoragecratio;

        for (int l=0; l!=2; ++l){
          soilcarr[0]->SOM[l+1] += carbondrawnLeaf*pft.leaflitterfrc[l];
          soilcarr[0]->SOM[l+1] += carbondrawnStem*pft.stemlitterfrc[l];
        }

        // root litter transfer
        for (int k=0; k!=ncells; ++k) {
          double carbondrawn = mort*pft.BRootSoil[k];
          soilcarr[k]->SOM[0] += carbondrawn*rootstoragecratio;
          //pft.Bstore = pft.Bstore - carbondrawn*rootstoragecratio;

          for (int l=0; l!=2; ++l){
            soilcarr[k]->SOM[l+1] = soilcarr[k]->SOM[l+1] + carbondrawn*pft.rootlitterfrc[l];
          }
          pft.BRootSoil[k] = pft.BRootSoil[k] * (1.0 - mort);
        }

        pft.Bleaf *= 1.0 - mort;
        pft.Broot *= 1.0 - mort;
        pft.Bstem *= 1.0 - mort;
        pft.Bstore *= 1.0 - mort;
        pft.totalBiomass = pft.Bleaf + pft.Broot + pft.Bstem + pft.Bstore;
      }

      //=============================================================================================
      // annual variable zeros
      if (std::fmod(t_days, 365.25) > std::fmod(t_days + dt_days, 365.25)) {
        // year rolled over

        // estimate the maximum leaf area index
        for (int i=0; i!=10; ++i){
          pft.annCBalance[i] = pft.annCBalance[i] - 1.0 / pft.SLA;
          pft.annCBalance[i] = pft.annCBalance[i] - 1.0 / pft.SLA
              * pft.root2leafrespratio * 1.0 / (pft.rootlongevity);
          pft.annCBalance[i] = pft.annCBalance[i] - 1.0 / pft.SLA
              * pft.stem2leafrespratio * 1.0 / (pft.stemlongevity);
        }
        pft.maxLAI = 0;

        for (int i=0; i!=10; ++i){
          if (pft.annCBalance[i] > 0.0) pft.maxLAI = i+1;
        }
        if (pft.maxLAI == 0) pft.maxLAI = 1;
        if (pft.pft_type == "moss") { 
            pft.maxLAI = 1; // for moss, no more than 2 leaf layers
        }

        // zero out annual variables
        for (int i=0; i!=max_leaf_layers; ++i){
          pft.annCBalance[i] = 0;
        }
        pft.annNPP = 0.0;
      }

      //===============================================================================
      // Check root mass balance
      pft.AssertRootBalance_or_die();

      //--------------------------------------------------------------------------------------
      // transfer to the liteter pool for dead vegetations, based on turn over rates
      // leaf turn over
      stemstoragecratio = pft.Bstore / ((pft.Bleaf + pft.Bleafmemory)*pft.storagecleaf2sw + pft.Bstem + pft.Broot*pft.storagecroot2sw);
      leafstoragecratio = pft.storagecleaf2sw*stemstoragecratio;
      double  rootstoragecratio = pft.storagecroot2sw*stemstoragecratio;

      // AFFECTED BY DT! --etc
      double turnoverLeaf = (!pft.evergreen) ? 0. : dt_days / (pft.leaflongevity * 365.25);
      double carbondrawnLeaf = turnoverLeaf*pft.Bleaf;
      AMANZI_ASSERT(turnoverLeaf <= 0.9);

      double turnoverStem = dt_days / (pft.stemlongevity * 365.25);
      double carbondrawnStem = turnoverStem*pft.Bstem;
      AMANZI_ASSERT(turnoverStem <= 0.9);

      double turnoverRoot = dt_days / (pft.rootlongevity * 365.25);
      AMANZI_ASSERT(turnoverRoot <= 0.9);

      double storagecdrawn = carbondrawnLeaf*leafstoragecratio*(1.0-pft.storagecRspFrc);
      soilcarr[0]->SOM[0] += storagecdrawn;
      pft.Bstore -= storagecdrawn;
      storagecdrawn = carbondrawnStem*stemstoragecratio*(1.0 - pft.storagecRspFrc);
      soilcarr[0]->SOM[0] += storagecdrawn;
      pft.Bstore -= storagecdrawn;

      for (int l=0; l!=2; ++l) {
          soilcarr[0]->SOM[l+1] += carbondrawnLeaf*pft.leaflitterfrc[l];
          soilcarr[0]->SOM[l+1] += carbondrawnStem*pft.stemlitterfrc[l];
      }

      //===============================================================================
      // Check root mass balance  --unnecessary.... it hasn't changed! --etc
      pft.AssertRootBalance_or_die();

      // root litter transfer
      for (int k=0; k!=ncells; ++k) {
          double carbondrawn = turnoverRoot *pft.BRootSoil[k];
          storagecdrawn = carbondrawn*rootstoragecratio*(1.0 - pft.storagecRspFrc);
          soilcarr[k]->SOM[0] += storagecdrawn;
          pft.Bstore -= storagecdrawn;
          for (int l=0; l!=2; ++l) {
            soilcarr[k]->SOM[l+1] += carbondrawn*pft.rootlitterfrc[l];
            AMANZI_ASSERT(pft.BRootSoil[k] >= 0.0);
          }

          pft.BRootSoil[k] = pft.BRootSoil[k] * (1.0 - turnoverRoot);
       }
       
      pft.Bleaf *= 1.0 - turnoverLeaf;
      pft.Broot *= 1.0 - turnoverRoot;
      pft.Bstem *= 1.0 - turnoverStem;
      pft.totalBiomass = pft.Bleaf + pft.Broot + pft.Bstem + pft.Bstore;

      //===============================================================================
      // Check root mass balance
      pft.AssertRootBalance_or_die();

      //=============================================================================
      //distribution transpirations
      for (int k=0; k!=ncells; ++k) {
        rootFrac = pft.Broot > 0. ? pft.BRootSoil[k] / pft.Broot : 0.;
        TransArr[k] += pft.ET * rootFrac/(gridarea * SoilThicknessArr[k] * dt);  //unit: kg H2O/m3/s
      }

    } else {  //biomass check

      if (std::fmod(t_days, 365.25) > std::fmod(t_days + dt_days, 365.25)) {
        // year rolled over
        // annual setup
        for (int i=0; i!=max_leaf_layers; ++i) {
          pft.annCBalance[i] = 0;
        }

        //seed rain
        pft.Bleaf = pft.seedrainlai*gridarea / pft.SLA;
        pft.Bstem = pft.Bleaf / pft.leaf2stemratio;
        pft.Broot = pft.Bleaf / pft.leaf2rootratio;
        if (!pft.evergreen) {
          pft.Bleafmemory = pft.Bleaf;
          pft.Bleaf = 0.0;
          pft.GDD = 0.0;
          pft.leafondaysi = 0.0;
          pft.leafstatus = 1;
        }
        pft.lai = pft.Bleaf*pft.SLA / gridarea;
        pft.laimemory = pft.Bleafmemory*pft.SLA / gridarea;
        pft.totalBiomass = pft.Bleaf + pft.Broot + pft.Bstem + pft.Bstore;
        pft.BRootSoil[0] = pft.Broot;
        double storagecleaf = (pft.Bleaf + pft.Bleafmemory)* pft.tar_leafstorageccon;
        pft.Bstore += storagecleaf;
        pft.Bstore += pft.Bstem / ((pft.Bleaf + pft.Bleafmemory)* pft.storagecleaf2sw )*storagecleaf;
        pft.Bstore += pft.Broot / (pft.Bleaf + pft.Bleafmemory)*storagecleaf * pft.storagecroot2sw / pft.storagecleaf2sw;
      } //end of year check 
    } //biomass check
  } // loop for different PFTs

  //---------------------------------------------------------------------------------
  //calculate shaded radiations for soil
  double radi = met.qSWin;
  for (std::vector<Teuchos::RCP<PFT> >::iterator pft_iter=pftarr.begin();
       pft_iter!=pftarr.end(); ++pft_iter) {
    //    std::cout << "wtf: (" << (*pft_iter)->pft_type << ") " << (*pft_iter)->LER << ", " << (*pft_iter)->lai << ", " << radi << std::endl;
    radi *= std::exp(-(*pft_iter)->LER * (*pft_iter)->lai);
   }
  sw_shaded = radi;

      

 //=========================================================================
  // do soil decomposition
  for (int k=0; k!=ncells; ++k) {
    double TFactor = TEffectsQ10(2.0, SoilTArr[k] - 273.15, 25.0);
    SoilCO2Arr[k] = 0.0;
    double WFactor;

    double soil_wp = std::max(std::min( (SoilWPArr[k] - p_atm) / 1.e6, wp_max), wp_min);
    if (soil_wp == wp_min) {
      WFactor = 0.0;
    } else {
      WFactor = std::log(wp_min / soil_wp) / std::log(wp_min / wp_max);
    }

    double DFactor = std::exp(-SoilDArr[k] / 0.5);
    int nPools = soilcarr[k]->params->nPools;
    std::vector<double> SOMConvt(nPools);

    for (int l=0; l!=nPools; ++l){
      // AFFECTED BY DT --etc
      double turnover = dt_days * WFactor*TFactor*DFactor / (soilcarr[k]->params->TurnoverRates[l] * 365.25);
      AMANZI_ASSERT(turnover <= 0.9);

      double SOMDecomp = soilcarr[k]->SOM[l] * turnover;
      soilcarr[k]->SOM[l] *= 1.0 - turnover;

      SoilCO2Arr[k] += SOMDecomp * soilcarr[k]->params->RespF[l];
      SOMConvt[l] = SOMDecomp * (1.0 - soilcarr[k]->params->RespF[l]);
    }

    double totalSoilCz = 0.0;
    for (int l=0; l!=nPools; ++l) {
      double totalConvrtC = 0.0;
      for (int m=0; m!=nPools; ++m) {
        if (m != l) {
          totalConvrtC += SOMConvt[m]* soilcarr[k]->params->Tij[m][l]; // *SOMConvt[l];
        }
      }

      soilcarr[k]->SOM[l] += totalConvrtC;
      totalSoilCz += soilcarr[k]->SOM[l];
    }
  }

  //================================================
  //do vertical diffusion
  Cryoturbate(dt_days, SoilTArr, SoilDArr, SoilThicknessArr, soilcarr, cryoturbation_coef);
  return;

  }

  // Cryoturbation -- move the carbon around via diffusion
  void Cryoturbate(double dt,
		   const Epetra_SerialDenseVector& SoilTArr,
		   const Epetra_SerialDenseVector& SoilDArr,
		   const Epetra_SerialDenseVector& SoilThicknessArr,
		   std::vector<Teuchos::RCP<SoilCarbon> >& soilcarr,
		   double diffusion_coef) {
    std::vector<double> diffusion_coefs(soilcarr[0]->nPools, diffusion_coef);
    Cryoturbate(dt, SoilTArr, SoilDArr, SoilThicknessArr, soilcarr, diffusion_coefs);
  }


  // Cryoturbation -- move the carbon around via diffusion
  void Cryoturbate(double dt,
		   const Epetra_SerialDenseVector& SoilTArr,
		   const Epetra_SerialDenseVector& SoilDArr,
		   const Epetra_SerialDenseVector& SoilThicknessArr,
		   std::vector<Teuchos::RCP<SoilCarbon> >& soilcarr,
		   std::vector<double>& diffusion_coefs) {
    int ncells = SoilTArr.Length();
    // only cryoturbate unfrozen soil
    int k_frozen = PermafrostDepthIndex(SoilTArr, 273.15);

    // fast and dirty diffusion
    int npools = soilcarr[0]->nPools;
    Epetra_SerialDenseVector dC_up(npools);
    Epetra_SerialDenseVector dC_dn(npools);
    std::vector<Epetra_SerialDenseVector> dC(k_frozen, Epetra_SerialDenseVector(npools));

    int k = 0;
    while (k < k_frozen) {
      Epetra_SerialDenseVector& C = soilcarr[k]->SOM;

      // dC/dz on the face above
      if (k == 0) {
	// dC_up initialized to zero already
      } else {
	Epetra_SerialDenseVector& C_up = soilcarr[k-1]->SOM;
	double dz_up = SoilDArr[k] - SoilDArr[k-1];
	for (int l=0; l!=npools; ++l) {
	  dC_up[l] = (C[l] - C_up[l]) / dz_up;
	}
      }

      // dC/dz on the face below
      if ((k == ncells-1) || k+1 >= k_frozen) {
	// boundary case
	for (int l=0; l!=npools; ++l) {
	  dC_dn[l] = 0;
	}
      } else {
	Epetra_SerialDenseVector& C_dn = soilcarr[k+1]->SOM;
	double dz_dn = SoilDArr[k+1] - SoilDArr[k];
	for (int l=0; l!=npools; ++l) {
	  dC_dn[l] = (C_dn[l] - C[l]) / dz_dn;
	}
      }

      // dC = dt * D * (dC/dz_below - dC/dz_above) / dz
      double dz = SoilThicknessArr[k];
      for (int l=0; l!=npools; ++l) {
	dC[k][l] = dt * diffusion_coefs[l] / dz * (dC_dn[l] - dC_up[l]);
      }

      // increment
      k++;
    }
  
    // Now that dC is calculated everywhere, update carbon pools
    k = 0;
    while (k < k_frozen) {
      Epetra_SerialDenseVector& C_k = soilcarr[k]->SOM;
      Epetra_SerialDenseVector& dC_k = dC[k];
    
      for (int l=0; l!=npools; ++l) {
	C_k[l] += dC_k[l];
      }
      k++;
    }
  }


} // namespace BGC
} // namespace Amanzi

