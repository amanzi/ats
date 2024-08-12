/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity of lake model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

ThermalConductivityEvaluator::ThermalConductivityEvaluator(
    Teuchos::ParameterList& plist) :
        SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key",
        "surface-thermal_conductivity");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  // -- depth
  depth_key_ = Keys::readKey(plist_, domain_name, "water depth", "water_depth");
  dependencies_.insert(depth_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  //  K_liq_ = sublist.get<double>("thermal conductivity of water [W/(m-K)]", 0.58);
  //  K_ice_ = sublist.get<double>("thermal conductivity of ice [W/(m-K)]", 2.18);
  //  min_K_ = sublist.get<double>("minimum thermal conductivity", 1.e-14);

  // later: read these parameters from xml
  K_max_    = 150.; // [W/(m * K)]
  K_0_      = 0.561;
  V_wind_   = 5.; // [m/s]
  V_wind_0_ = 40.;
}


ThermalConductivityEvaluator::ThermalConductivityEvaluator(
    const ThermalConductivityEvaluator& other) :
        SecondaryVariableFieldEvaluator(other),
        K_max_(other.K_max_),
        K_0_(other.K_0_),
        V_wind_(other.V_wind_),
        V_wind_0_(other.V_wind_0_),
        temperature_key_(other.temperature_key_),
        depth_key_(other.depth_key_){}
//    uf_key_(other.uf_key_),
//    height_key_(other.height_key_),
//    K_liq_(other.K_liq_),
//    K_ice_(other.K_ice_),
//    min_K_(other.min_K_) {}


Teuchos::RCP<FieldEvaluator>
ThermalConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityEvaluator(*this));
}

void ThermalConductivityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  ice_cover_ = false; // first always assume that there is no ice

  double z_ice = 1.0;
  double z_w   = 1.0;

  double lambda_ice = 2.2;
  double lambda_w   = 0.561; //3.*0.561; //1.5
  // lambda_ice = lambda_w;

  // read parameters from the met data
  Teuchos::ParameterList& param_list = plist_.sublist("parameters");
  Amanzi::FunctionFactory fac;
  Teuchos::RCP<Amanzi::Function> uwind_func_ = Teuchos::rcp(fac.Create(param_list.sublist("wind u-speed")));
  Teuchos::RCP<Amanzi::Function> vwind_func_ = Teuchos::rcp(fac.Create(param_list.sublist("wind v-speed")));

  std::vector<double> args(1);
  args[0] = S->time();
  double uwind = (*uwind_func_)(args);
  double vwind = (*vwind_func_)(args);
  V_wind_ = sqrt(uwind*uwind+vwind*vwind);

  // std::cout << "V_wind = " << V_wind_ << std::endl;

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  // get depth
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  /// WHY ISN'T IT INITIALIZED ????

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    //      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp,false);
    //      const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
    const Epetra_MultiVector& h = *depth->ViewComponent(*comp,false);

    Key domain = Keys::getDomain(my_key_);

    // Set up my dependencies.
    std::string domain_name = Keys::getDomain(my_key_);

    const Epetra_MultiVector& cv =
      *S->GetFieldData(Keys::getKey(domain_name,"cell_volume"))->ViewComponent("cell",false);

    int ncomp = result->size(*comp, false);

    //      std::cout << "ncomp in lambda = " << ncomp << std::endl;

    int i_ice_max = 0;
    int i_ice_min = ncomp;

    for (int i=0; i!=ncomp; ++i) {
      if (temp_v[0][i] < 273.15) { // check if there is ice cover
        ice_cover_ = true;
        i_ice_max = i;
      }
    } // i

    for (int i=ncomp-1; i!=-1; --i) {
      if (temp_v[0][i] < 273.15) { // check if there is ice cover
        ice_cover_ = true;
        i_ice_min = i;
      }
    } // i

    int d_thawed = ncomp-1-i_ice_max;  // thickness of thawed layer [cells]
    int d_ice = i_ice_max-i_ice_min+1; // thickness of ice layer [cells]

    if (ice_cover_) {
      const AmanziGeometry::Point& zci = mesh->cell_centroid(i_ice_max);
      z_ice = zci[2];

      const AmanziGeometry::Point& zcw = mesh->cell_centroid(i_ice_min);
      z_w = zcw[2];
    }

    for (int i=0; i!=ncomp; ++i) {
      if (temp_v[0][i] < 273.15) { // this cell is in ice layer
        result_v[0][i] = lambda_ice; //0.5*(lambda_ice + result_v[0][i]);
      }
      else { // this cell is in water layer
        double lambda_tmp = result_v[0][i];
        if (ice_cover_) {
          lambda_tmp = lambda_w; //1.5; 
        } else {
          lambda_tmp = lambda_w; //10.*K_0_ + V_wind_/V_wind_0_*(K_max_ - K_0_);
        }
        // std::cout << "lambda_tmp = " << lambda_tmp << std::endl;
        // lambda_tmp = 1.5;
        result_v[0][i] = lambda_tmp; //0.5*(lambda_tmp + result_v[0][i]);
      }
    } // i

    std::vector<double> lambda_new(ncomp);

    for (int i=0; i!=ncomp; ++i) {
      lambda_new[i] = result_v[0][i];
    }

    if (ice_cover_ && d_thawed > 0) {

//      std::cout << "i_ice_min = " << i_ice_min << std::endl;
//      std::cout << "i_ice_max = " << i_ice_max << std::endl;
//      std::cout << "d_thawed = " << d_thawed << std::endl;
//      std::cout << "d_ice    = " << d_ice << std::endl;

//      std::cout << "lambda before swap " << std::endl;
      for (int i=ncomp-1; i!=-1; --i) {
//        std::cout << "result_v[0][" << i << "] = " << result_v[0][i] << std::endl;
      }

      // if thawing occured at the top, swap cells
      for (int i=0; i < d_ice; ++i) { // push ice to the surface
//        std::cout << "copy cell " << i_ice_max-i << " to " << ncomp-1-i << std::endl;
        lambda_new[ncomp-1-i] = result_v[0][i_ice_max-i];
      }
      for (int i=0; i < d_thawed; ++i) { // push water to the bottom
        lambda_new[i_ice_min+i] = result_v[0][i_ice_max+1+i];
//        std::cout << "copy cell " << i_ice_max+1+i << " to " << i_ice_min+i << std::endl;
      } // i

      for (int i=0; i!=ncomp; ++i) {
      //  result_v[0][i] = lambda_new[i];
      }

//      std::cout << "lambda after swap " << std::endl;
      for (int i=ncomp-1; i!=-1; --i) {
//        std::cout << "result_v[0][" << i << "] = " << result_v[0][i] << std::endl;
        if (result_v[0][i] < 0.) exit(0);
      }

//      exit(0);

    }

    /*
    std::vector<double> temp_new(ncomp); // new temperatures for swapping the cells

    for (int i=0; i!=ncomp; ++i) {
      temp_new[i] = temp_v[0][i];
    }

    if (ice_cover_ && d_thawed > 0) {

      std::cout << "i_ice_min = " << i_ice_min << std::endl;
      std::cout << "i_ice_max = " << i_ice_max << std::endl;
      std::cout << "d_thawed = " << d_thawed << std::endl;
      std::cout << "d_ice    = " << d_ice << std::endl;


      std::cout << "Temperature before swap " << std::endl;
      for (int i=ncomp-1; i!=-1; --i) {
        std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
      }

      // if thawing occured at the top, swap cells
      for (int i=0; i < d_ice; ++i) { // push ice to the surface
        std::cout << "copy cell " << i_ice_max-i << " to " << ncomp-1-i << std::endl;
        temp_new[ncomp-1-i] = temp_v[0][i_ice_max-i];
      }
      for (int i=0; i < d_thawed; ++i) { // push water to the bottom
        temp_new[i_ice_min+i] = temp_v[0][i_ice_max+1+i];
        std::cout << "copy cell " << i_ice_max+1+i << " to " << i << std::endl;
      } // i

      for (int i=0; i!=ncomp; ++i) {
        temp_v[0][i] = temp_new[i];
      }

      std::cout << "Temperature after swap " << std::endl;
      for (int i=ncomp-1; i!=-1; --i) {
        std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
      }

  //    if (d_ice < d_thawed) exit(0);

    }

    */

    // simple interpolation between neighboring cells
    std::vector<double> lambda(ncomp);

    for (int i=0; i!=ncomp; ++i) {
      lambda[i] = result_v[0][i];
    } // i

    // lambda[0] = result_v[0][0];
    // for (int i=1; i!=ncomp; ++i) {
    //   const AmanziGeometry::Point& zc = mesh->cell_centroid(i);
    //   lambda[i] = 0.5*(result_v[0][i-1]+result_v[0][i]);
    // } // i

  //  lambda[0] = result_v[0][0];
  //  lambda[ncomp-1] = result_v[0][ncomp-1];
  //  for (int i=1; i!=ncomp-1; ++i) {
  //    const AmanziGeometry::Point& zc = mesh->cell_centroid(i);
  //    lambda[i] = 1./3.*(result_v[0][i-1]+result_v[0][i]+result_v[0][i+1]);
  //  } // i

    // std::cout << "depth in eval = " << h[0][0] << ", squared = " << h[0][0]*h[0][0] << std::endl;

    for (int i=0; i!=ncomp; ++i)  {
      // std::cout << "lambda[i] = " << lambda[i] << std::endl;
      double hh = h[0][0]; //(h[0][0] == 0.) ? 1.5 : h[0][0]; // FoxDen
      // double hh = (h[0][0] == 0.) ? 3.0 : h[0][0]; // Atqasuk
      // double hh = (h[0][0] == 0.) ? 24.0 : h[0][0]; // Toolik
      // std::cout << "hh = " << hh << std::endl;
      result_v[0][i] = lambda[i]/(hh*hh); ///cv[0][0];
      // std::cout << "result = " << result_v[0][i] << std::endl;
    }

  }

}


void ThermalConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
