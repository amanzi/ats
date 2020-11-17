/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov)
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "dbc.hh"
#include "pc_ice_water.hh"

namespace Amanzi {
namespace Flow {


/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
PCIceWater::PCIceWater(Teuchos::ParameterList& pc_plist) :
    pc_plist_(pc_plist) {
  InitializeFromPlist_();
};

double PCIceWater::CapillaryPressure(double T, double dens) {
  double pc;
  if (halfwidth_ > 0) {
    double dx = halfwidth_;
    double alpha = gamma_*dens/T0_;
    double a = alpha* dx/4.;
    double b = -alpha/2.;
    double c = alpha/(4.*dx);
    double x = T-T0_;

    if(x < -dx) {
      pc = -alpha* x;
    } else if( x > dx ) {
      pc = 0.;
    } else {
      pc = a + b*x + c*x*x;
    }
  } else {
    pc = T < T0_ ? gamma_ * dens * (T0_ - T)/T0_ : 0.;
  }
  return pc;
};

double PCIceWater::DCapillaryPressureDT(double T, double dens) {
  double dpc;
  if (halfwidth_ > 0.) {
    double dx = halfwidth_;
    double alpha = gamma_*dens/T0_;
    double a = alpha* dx /4.;
    double b = -1.0*alpha/2.;
    double c = alpha/(4.*dx);
    double x = T-T0_;

    if(x < -dx) {
      dpc = -alpha ;
    } else if( x > dx ) {
      dpc = 0.;
    } else {
      dpc =  b + 2.*c*x;
    }
  } else {
    dpc = T < T0_ ? -gamma_ * dens / T0_ : 0.;
  }
  return dpc;
};

double PCIceWater::DCapillaryPressureDRho(double T, double dens) {
  double dpc;
  if (halfwidth_ > 0.) {
    double dx = halfwidth_;
    double dalpha = gamma_/T0_;
    double a = dalpha* dx/4.;
    double b = -1.0*dalpha/2.;
    double c = dalpha/(4.*dx);
    double x = T-T0_;

    if(x < -dx) {
      dpc = -dalpha * x ;
    } else if( x > dx ) {
      dpc = 0.;
    } else {
      dpc =  a + b *x + c* x*x;
    }
  } else {
    dpc = T < T0_ ? gamma_ * (T0_ - T)/T0_ : 0.;
  }
  return dpc;
};


void PCIceWater::InitializeFromPlist_() {
  sigma_ice_liq_ = pc_plist_.get<double>("interfacial tension ice-water [mN m^-1]", 33.1);
  sigma_gas_liq_ = pc_plist_.get<double>("interfacial tension air-water [mN m^-1]", 72.7);
  T0_ = pc_plist_.get<double>("reference temperature [K]", 273.15);
  halfwidth_  = pc_plist_.get<double>("smoothing width [K]", -1.)/2.;

  if (pc_plist_.isParameter("latent heat [J mol^-1]")) {
    heat_fusion_ = pc_plist_.get<double>("latent heat [J mol^-1]");
    molar_basis_ = true;
  } else {
    heat_fusion_ = pc_plist_.get<double>("latent heat [J kg^-1]", 3.34e5);
    molar_basis_ = false;
  }

  gamma_ = sigma_gas_liq_/sigma_ice_liq_ * heat_fusion_;
};

} // namespace
} // namespace
