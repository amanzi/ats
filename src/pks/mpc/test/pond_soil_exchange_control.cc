#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <limits>
#include "pond_soil_exchange_control.hh"
int main(){
 using Amanzi::PondSoilExchangeLimit;
 auto f=[](double h,double cp,double sensible,double cs,double g){return PondSoilExchangeLimit(h,cp,sensible,cs,g,1e-5,10.,4.,3600.);};
 const auto x=f(.000616,.993487,.002461,19.411,9.1e-6);
 assert(x.contrast>19. && x.contrast<20.);
 assert(x.tau>260. && x.tau<280. && x.apparent_tau>1e5);
 assert(x.timestep==900.);
 assert(std::isinf(f(0.,.993487,.002461,19.411,9.1e-6).timestep));
 assert(std::isinf(f(.000616,.993487,.002461,.993487,9.1e-6).timestep));
 assert(std::isinf(f(.000616,.993487,.002461,19.411,0.).timestep));
 assert(std::isinf(f(.000616,1.003487,.012461,19.411,9.1e-6).timestep));
 const auto y=f(.000616,993487.,2461.,19411000.,9.1);
 assert(std::abs(x.tau-y.tau)<1e-10 && x.timestep==y.timestep);
 assert(f(.000616,.993487,.002461,19.411,9.1e-7).timestep>=x.timestep);
 assert(std::isinf(f(.000616,.993487,.002461,19.411,std::numeric_limits<double>::infinity()).timestep));
 assert(std::isinf(f(.000616,std::numeric_limits<double>::quiet_NaN(),.002461,19.411,9.1e-6).timestep));
 assert(std::isinf(f(.000616,-1.,.002461,19.411,9.1e-6).timestep));
}
