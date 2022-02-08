/*
  This is one of the flow component of the ATS physics code.
  License: BSD
  Authors:
  F.-M. Yuan (yuanf@ornl.gov)
*/

#include <cmath>
#include <fstream>
#include "dbc.hh"
#include "errors.hh"
#include "Spline.hh"

#include "wrm_clapp_hornberger.hh"

namespace Amanzi {
namespace Flow {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMClappHornberger::WRMClappHornberger(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};


/* ******************************************************************
* Relative permeability formula: input is liquid saturation.
* Krel = Se**(3+2*bsw)
*
****************************************************************** */
double WRMClappHornberger::k_relative(double s) {
  double se = (s - sr_)/(1.0-sr_);
  double sex = std::max(0.0, (sx_ - sr_)/(1.0-sr_));
  double se0 = std::max(0.0, (s0_ - sr_)/(1.0-sr_));
  if (se <= sex) {
    return 0.0;
  } else if (se >= 1.0) {
    return 1.0;
  } else if (se > se0) {
    return fit_kr0_(s);

  } else {
    return std::pow(se, (3.0+2.0*bsw_));
  }
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double WRMClappHornberger::d_k_relative(double s) {
  double se = (s - sr_)/(1.0-sr_);
  double sex = std::max(0.0, (sx_ - sr_)/(1.0-sr_));
  double se0 = std::max(0.0, (s0_ - sr_)/(1.0-sr_));

  if (se <=sex || se>=1.0) {
    return 0.0;
  } else if (se > se0) {
    return fit_kr0_.Derivative(s);

  } else {
    double dk = (3.0+2.0*bsw_)*std::pow(se, (2.0+2.0*bsw_));
    return dk/(1.0-sr_);

  }
}


/* ******************************************************************
 * Saturation formula: inverse of pc = smpsat * S**(-bsw)
 ****************************************************************** */
double WRMClappHornberger::saturation(double pc) {

  double pr = std::max(pc/smpsat_, 1.0);
  double se = std::pow(pr, -1.0/bsw_);

  if (pc < pc0_) {

    if (pc <= 0.0 || (pc0_ <= smpsat_ || s0_ >= 1.0)) {
      se = 1.0;
      return se*(1.0-sr_)+sr_;

    } else {
      // near-saturation (inflection point Si ~ 1.0) parabolic decreasing to zero
      if (!wet_step_smoothing_) {
        // pc = -m*(se-n)*(se-1.0)=-m*(se^2-(n+1)*se+n)
        // se^2-(n+1)*se+(n+pc/m) = 0
        pr = (n_+1.0)*(n_+1.0)/4.0-(n_+pc/m_);
        if (pr<0.0) {  // checking only, mathematically not possible
    	    Errors::Message message("WRM: Clapp Hornberger Inflection point NOT right! ");
            std::cout<< "PSIi: "<<pc0_ << "Si: "<<s0_
        		  <<"m: "<<m_ <<"n: "<<n_ <<std::endl;
    	    Exceptions::amanzi_throw(message);
        } else {
          double se = std::sqrt(pr)+(n_+1.0)/2.0;
          return se*(1.0-sr_)+sr_;
        }
      }

    }
  }

  // inverse of sigmoid-smoothed PSI(sat) step function
  //if (wet_step_smoothing_){
  //  double s = se*(1.0-sr_)+sr_;
  //  double pc_factor = capillaryPressure(s)/pc;
  // double smoothing_factor = std::log(1.0/pc_factor-1.0)/(-sigmoid_k0_)+sigmoid_half0_;
  // se = smoothing_factor;
  //}

  return se * (1.0 - sr_) + sr_;

}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double WRMClappHornberger::d_saturation(double pc) {

  double pr = std::max(pc/smpsat_, 1.0);
  double dse = (std::pow(pr, -1.0/bsw_-1.0)) * (-1.0/bsw_);

  if (pc < pc0_) {

    if (pc <= 0.0 || (pc0_ <= smpsat_ || s0_ >= 1.0)) {
      return 0.0;

    } else {
      // near-saturation (inflection point PCi ~ 0.0) smoothing
      //
      if (!wet_step_smoothing_) {
        // pc = -m*(se-n)*(se-1.0)=-m*(se^2-(n+1)*se+n)
        // se^2-(n+1)*se+(n+pc/m) = 0
        double pr = (n_+1.0)*(n_+1.0)/4.0-(n_+pc/m_);
        if (pr<0.0) { // checking only, mathematically not possible
      	  Errors::Message message("WRM: Clapp Hornberger Inflection point NOT right! ");
          std::cout<< "PSIi: "<<pc0_ << "Si: "<<s0_
      	              <<"m: "<<m_ <<"n: "<<n_ <<std::endl;
      	  Exceptions::amanzi_throw(message);

        } else {
          double dse = 0.5*std::pow(pr, -0.5)*(-1.0/m_);
          return dse * (1.0 - sr_);
        }
      }
	}
  }

  return dse * (1.0 - sr_);

}

/* ******************************************************************
 * Pressure as a function of saturation: pc= smpsat * S**(-bsw)
 ****************************************************************** */
double WRMClappHornberger::capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 0.0);
  double pc;

  //s  ranges: -- 0 ---- sr ---- sx_ -- s0_ -- (1.0-sr_) -- 1.0 --
  //se       : --------- 0. -------------------- 1.0     ---------
  //pc       : -------- Inf ---- pcx_ --pc0_ --- 0.0     ---------
  //krel     : ----------------- 0.0  ---------- 1.0     ---------

  pc = smpsat_*(std::pow(se, -bsw_));

  if (s <= sx_ ) { // dry-end, if any
    if (sx_>sr_ && pcx_>0.0) {
      pc = pcx_;
    } else{
      pc = std::numeric_limits<double>::infinity();
    }

  } else if (s >= 1.0-sr_) { // over-saturated
    pc = 0.0;

  } else if (s > s0_) { // near-saturated: s0_ -- (1.0-sr_) --
	if (pc0_ <= smpsat_ || s0_ >= 1.0) {
      // neither smoothing nor inflection point decreasing
      pc = smpsat_;
    } else {
      if (!wet_step_smoothing_) {
        // near-saturation (inflection point Si ~ 1.0) parabolic decreasing to zero
        // pc = -m*(se-n)*(se-1.0)=-m*(se^2-(n+1)*se+n)
        pc = -m_*(se-n_)*(se-1.0);
      }
    }
  }

  // multiplicative of a sigmoid-smoothing step-function for all sat range
  if (wet_step_smoothing_) {
    double s_u =  se - sigmoid_half0_;
    double smoothing_factor = 1.0/(1.0+std::exp(-sigmoid_k0_*s_u));
    pc = pc*smoothing_factor;
  }

  return pc;
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double WRMClappHornberger::d_capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 0.0);
  double dpc;
  dpc = smpsat_*(-bsw_)*(std::pow(se, -bsw_-1.0))/(1.0-sr_);

  if ((sx_>sr_ && pcx_>0.0) && s <= sx_) {  // dry-end
    dpc = 0.0;

  } else if (s >= (1.0 - sr_)){  // over-saturated
    dpc = 0.0;

  } else if (s > s0_){ // near-saturated: (exclusive)s0_ -- (1.0-sr_) --
	if (pc0_ <= smpsat_ || s0_ >= 1.0) {
      dpc = 0.0;
	} else {
      //
      if (!wet_step_smoothing_) {
        // near-saturation (inflection point Si ~ 1.0) parabolic decreasing to zero
        // pc = -m*(se-n)*(se-1.0)=-m*(se^2-(n+1)*se+n)
        dpc = -m_*(2.0*se-n_-1.0)/(1.0-sr_);
      }
    }
  }

  // multiplicative of a sigmoid-smoothing step-function for all sat range
  if (wet_step_smoothing_) {
    double s_u =  se - sigmoid_half0_;
    double smoothing_factor = 1.0/(1.0+std::exp(-sigmoid_k0_*s_u));
    double d_smoothing = smoothing_factor * (1.0-smoothing_factor);
    wet_step_smoothing_ = false;   // reset it to false so that 'pc' below is the original
    double pc = capillaryPressure(s);
    dpc = dpc*smoothing_factor + pc*d_smoothing;
    wet_step_smoothing_ = true;    // set it back to true
  }

  return dpc;

}


void WRMClappHornberger::InitializeFromPlist_() {
  std::string fname = plist_.get<std::string>("Krel function name", "");
  if (fname != std::string("")) {
	Errors::Message message("WRM: Clapp Hornberger \"Krel function name\" CANNOT be user-defined");
	Exceptions::amanzi_throw(message);
  }

  smpsat_ = plist_.get<double>("Clapp Hornberger smpsat [Pa]", 4689.0); // Pa <-- 9.80665*478 mmH2O Head (loam soil, Clapp & Hornberger 1978)
  if (smpsat_ <= 0.0) {
	Errors::Message message("WRM: Clapp Hornberger \"Clapp Hornberger smpsat [Pa]\" CANNOT be less than 0.0");
	Exceptions::amanzi_throw(message);
  }
  bsw_ = plist_.get<double>("Clapp Hornberger bsw [-]", 5.39);          // (loam soil, Clapp & Hornberger 1978)
  if (bsw_ <= 1.0) {
	Errors::Message message("WRM: Clapp Hornberger \"Clapp Hornberger bsw [-]\" CANNOT be less than 1.0");
	Exceptions::amanzi_throw(message);
  }
  sr_ = plist_.get<double>("residual saturation [-]", 0.0);
  
  // Near-saturation smoothing, i.e. from 1.0 ~ Wi (inflection point, 0.92 in Clapp & Hornberger 1978)
  s0_ = 1.0 - plist_.get<double>("near-saturation inflection point interval [saturation]", 0.08);  // if user-input as non-positive, it will be off
  pc0_ = 0.0;
  wet_step_smoothing_ = false;
  if (s0_ < 1.0) {
    pc0_ = capillaryPressure(s0_);

  } else {

	// linear derivative-smoothing approach
    pc0_ = plist_.get<double>("near-saturation smoothing starting point [Pa]", 0.0);   // if user-input as less than 'smpsat', it will be off
    if (pc0_ > smpsat_) {
      s0_ = saturation(pc0_);

      double H1 = 1.0;
      double H0 = 0.0;
      double se0 = (s0_-sr_)/(1.0-sr_);
      sigmoid_one_eps_ = 1.e-4;
      sigmoid_half0_ = (se0+1.0)/2.0;            // sigmoid x ranges from eff. Sat of se0 to 1.0
      sigmoid_k0_ = std::log((H0+sigmoid_one_eps_)/(H1-sigmoid_one_eps_))
                    *2.0/(1.0-se0);             //

      //wet_step_smoothing_ = true; // NOT YET ready
    }
  }
  // clapp-hornberger (1978) approach
  if (s0_ < 1.0 && pc0_ > smpsat_) {
    m_ = pc0_/(1.0-s0_)/(1.0-s0_)-pc0_*bsw_/s0_/(1.0-s0_);
    n_ = 2.0*s0_-pc0_*bsw_/m_/s0_-1.0;

    // a simple checking on parameters chosen by user
    double se_sat = 0.9999;
    double pc_sat = -m_*(se_sat-n_)*(se_sat-1.0);
    if (pc_sat<0.0) {
      Errors::Message message("WRM: For the Clapp Hornberger Eqs, near-saturation inflection point NOT good!"
    		  " - probably too far away from saturaton");
  	  Exceptions::amanzi_throw(message);

    }

    // smoothing kr as well
    fit_kr0_.Setup(s0_, k_relative(s0_), d_k_relative(s0_),
  	       1.0, 1.0, 0.0);

  }


  // IF on, dry-end trucation or smoothing (TODO)
  pcx_ = plist_.get<double>("dry-end smoothing starting point [Pa]", 1.0e10); // user-input as non-positive, it will be off
  dry_step_smoothing_ = false;
  if (pcx_<=0.0) {
    sx_ = sr_;
  } else {
    sx_ = saturation(pcx_);  // mathematically >= sr_
  }

  // the following is for testing if wrm models coded correctly
  // by writint out a ascii csv file for processing and inspecting
//#define TEST_WRM_ClappHornberger
#ifdef TEST_WRM_ClappHornberger
    // generate a few text files for checking water retention curves
    // note: if there are more than 1 materials, this will be the last one
    std::ofstream txtout("wrm_clapphornberger_check.txt");
    std::streambuf *coutbuf = std::cout.rdbuf();  //before redirect, save old cout
    std::cout.rdbuf(txtout.rdbuf());              //redirect std::cout to txtout

    int count = 1000;
    double sat;
    double pc;
    double krel;

    double eps=1e-8;
    bool warned = false;
    std::cout << "Clapp-Hornberger WRM " << std::endl;
    std::cout << "smpsat: "<< smpsat_<<std::endl;
    std::cout << "bsw: "<< bsw_<<std::endl;
    std::cout << "Sr: "<< sr_<<std::endl;
    std::cout << "inflection point: Si "<< s0_<<std::endl;
    std::cout << "inflection point: PSIi "<< pc0_<<std::endl;
    if(wet_step_smoothing_) {
    	std::cout << "wet_step_smoothing sigmoid-k0: "<< sigmoid_k0_<<std::endl;
    }
    std::cout << "dry-end starting point: Sx "<< sx_<<std::endl;
    std::cout << "dry-end starting point: PSIx "<< pcx_<<std::endl;
    std::cout <<std::endl <<"-- i ---- Sat ---- PSI(Pa) ---- dPSI/dSat ---- Krel ---- dKrel ---- "<<std::endl;

    for (int i=0; i!=count; ++i) {
      sat = i/(count-1.0)*(1.0-sx_)+sx_;
      pc = capillaryPressure(sat);
      //pc  = std::pow(10, i/100.0);
      //sat = saturation(pc);

      if (abs(sat - saturation(pc)) > eps) {
        std::cout << "ERROR: sat != wrm.saturation(wrm.capillaryPressure(sat)) "
      		<<sat<< "- "<< saturation(pc)<< std::endl;
        warned = true;

      } else if (abs(pc - capillaryPressure(sat)) > eps*1.e5) {
        std::cout << "ERROR: pc[kPa] != wrm.capillaryPressure(wrm.saturation(pc))[kPa] "
        		<<pc*1.0e-3<< "- "<< capillaryPressure(sat)*1.0e-3<< std::endl;
        warned = true;
      }
      krel = k_relative(sat);

      std::cout << i <<", " << sat << ", " <<pc <<", "<< d_capillaryPressure(sat)
    		  <<", "<<krel << ", "<< d_k_relative(sat)<<std::endl;

    }
    std::cout<<std::endl;
    std::cout.rdbuf(coutbuf); //reset to standard output again

    if (warned) {
      Errors::Message message("WRM: Clapp Hornberger Eqs NOT right, please checking file: wrm_clapphornberger_check.txt ");
	  Exceptions::amanzi_throw(message);
	}

#endif


};

/* ******************************************************************
* Suction formula: input is liquid saturation.
****************************************************************** */
double WRMClappHornberger::suction_head(double s) {
  double pc = capillaryPressure(s);
  return -pc;
}


/* ******************************************************************
 * D suction_head / D saturation
 ****************************************************************** */
double WRMClappHornberger::d_suction_head(double s) {
  double dpc = d_capillaryPressure(s);
  return -dpc;
}

}  // namespace
}  // namespace
 
