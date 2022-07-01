/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS driver's parameter list method class for ATS-ELM interface

License: see $ATS_DIR/COPYRIGHT
Authors: Ethan Coon, Joe Beisman, F.-M. Yuan @ ORNL

Implementation (driver) of parameter list changing in the ats_elm_interface.
This is in contrast to directly over-ride State evaluator's values.
And must be called prior to model setup, but after read-in input *.xml

------------------------------------------------------------------------- */
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>

#include "Units.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "errors.hh"
#include "exceptions.hh"

// -----------------------------------------------------------------------------
//
// ats_elm interface parameter_list methods
//
// -----------------------------------------------------------------------------

#include "elm_ats_plist.hh"

namespace ATS {

elm_ats_plist::elm_ats_plist(Teuchos::RCP<Teuchos::ParameterList> &plist){
  //
  parameter_list_  = plist;
  elm_drv_plist_   = Teuchos::sublist(plist, "cycle driver");
  elm_pks_plist_   = Teuchos::sublist(plist, "PKs");
  elm_state_plist_ = Teuchos::sublist(plist, "state");

  subpk_plist_ = Teuchos::sublist(elm_pks_plist_, "flow");
  srfpk_plist_ = Teuchos::sublist(elm_pks_plist_, "overland flow");
}

// passing some of ELM data by resetting 'parameter list' in *.xml
// but have to do this prior to model setup
void elm_ats_plist::set_plist(elm_data &elmdata, const double start_ts, const double dt) {

  // (1) ELM surface grids/soil columns passing into ats mesh parameter lists
  //     if mesh type is 'generate mesh' from box-range
  plist_general_mesh_reset(elmdata);

  // (2) 'cycle driver' or 'coordinator' parameter-list
  plist_cycle_driver_reset(start_ts, dt);

  // (3) materials reset via plist
  plist_materials_reset(elmdata);

}

// ELM domain, surface-grids/soil profile, passing into ats via 'plist'
void elm_ats_plist::plist_general_mesh_reset(elm_data &elmdata, const bool elm_matched) {
  //
  int NX = elmdata.NX_;
  int NY = elmdata.NY_;
  int NZ = elmdata.NZ_;

  /*---------------------------------------------------------------------------------*/
  // (1) mesh, including 'surface' and 'domain'

  Teuchos::RCP<Teuchos::ParameterList> mesh_plist_ = Teuchos::sublist(parameter_list_, "mesh");
  Teuchos::RCP<Teuchos::ParameterList> domain_plist_ = Teuchos::sublist(mesh_plist_, "domain");

  if (domain_plist_->get<Teuchos::string>("mesh type") == "generate mesh") {

	  Teuchos::RCP<Teuchos::ParameterList>coords_plist_ = Teuchos::sublist(domain_plist_, "generate mesh parameters");

	  Teuchos::Array<double> coords_low;
	  Teuchos::Array<double> coords_high;
	  Teuchos::Array<double> coords_point;

	  //surface grids over-ride (TODO)

	  // modify vertical range
	  coords_low = coords_plist_->get<Teuchos::Array<double>>("domain low coordinate");
	  coords_high = coords_plist_->get<Teuchos::Array<double>>("domain high coordinate");

      // let's ats determine mesh, except for box size (ranges)
      coords_low[coords_low.size()-1] = elmdata.cols_verticesZ_[NZ-1];
      coords_high[coords_high.size()-1] = elmdata.cols_verticesZ_[0];
      Teuchos::Array<int> ncells = coords_plist_->get<Teuchos::Array<int>>("number of cells");
      ncells[0] = NX-1;
      ncells[1] = NY-1;
      ncells[2] = NZ-1;
      coords_plist_->set("number of cells", ncells);

      if (elm_matched) {
		//coords = elm_col_nodes;  //incorrect and not here
	  }

      coords_plist_->set("domain low coordinate", coords_low, "m");
	  coords_plist_->set("domain high coordinate", coords_high, "m");

	  /*---------------------------------------------------------------------------------*/
	  // (2) regions, including 'computational domain', 'surface domain', 'surface', 'bottom face', etc.
	  Teuchos::RCP<Teuchos::ParameterList> region_plist_ = Teuchos::sublist(parameter_list_, "regions");

	  Teuchos::RCP<Teuchos::ParameterList> compu_domain_ = Teuchos::sublist(region_plist_, "computational domain");
	  coords_plist_ = Teuchos::sublist(compu_domain_, "region: box");
	  coords_plist_->set("low coordinate", coords_low, "m");    // exactly same as 'mesh->domain'
	  coords_plist_->set("high coordinate", coords_high, "m");  // exactly same as 'mesh->domain'

	  Teuchos::RCP<Teuchos::ParameterList> sideset_surface_ = Teuchos::sublist(region_plist_, "surface");
	  coords_plist_ = Teuchos::sublist(sideset_surface_, "region: plane");
	  coords_point = coords_plist_->get<Teuchos::Array<double>>("point");
	  //coords_point[coords_point.size()-1] = elm_surf_gridsZ[0]; //TODO - this is what really needed
	  coords_point[coords_point.size()-1] = elmdata.cols_verticesZ_[0];
	  coords_plist_->set("point", coords_point, "m");

  }; // 'mesh type' is 'generate mesh'

}

// override whatever read-in from *.xml
void elm_ats_plist::plist_materials_reset(elm_data &elmdata) {

  // flow PK initial condition by water table
  Teuchos::RCP<Teuchos::ParameterList> pk_flow_initial_plist_ = Teuchos::sublist(
                    Teuchos::sublist(elm_pks_plist_, "flow"),  "initial condition");
    if (pk_flow_initial_plist_->get<double>("hydrostatic head [m]",-999.9)!=-999.9){
        double zwt = pk_flow_initial_plist_->get<double>("hydrostatic head [m]");
        zwt = elmdata.zwt_[0];   // (TODO) need to figure out how to pass multiple columns' data
        pk_flow_initial_plist_->set("hydrostatic head [m]", zwt);
    }
    
  // porosity, viscosity, & permibility in 'state -> evaluators'
  Teuchos::RCP<Teuchos::ParameterList> evals_plist_ = Teuchos::sublist(elm_state_plist_, "evaluators");
  Teuchos::RCP<Teuchos::ParameterList> state_constants_plist_ = Teuchos::sublist(elm_state_plist_, "initial conditions");

  // (1a) porosity
  auto base_poro_plist_ = Teuchos::sublist(
		  Teuchos::sublist(evals_plist_, "base_porosity"), "function");
  int c=0;
  for (auto& entry : *base_poro_plist_) {
    std::string layer_name = entry.first;

    Teuchos::RCP<Teuchos::ParameterList> layer_plist_ = Teuchos::sublist(
		                  Teuchos::sublist(
		                    Teuchos::sublist(base_poro_plist_, layer_name),"function"),
					    "function-constant");
    double poro = layer_plist_->get<double>("value");
    poro = elmdata.porosity_[c];
    layer_plist_->set("value", poro);
    c++;
  }

  // (1b) permeability
  auto visco_plist_ = Teuchos::sublist(evals_plist_, "viscosity_liquid");
  double visco = visco_plist_->get<double>("value");

  auto den_mass_plist_ = Teuchos::sublist(evals_plist_, "mass_density_liquid");
  double den_mass = den_mass_plist_->get<double>("value");

  auto gravity_plist_ = Teuchos::sublist(state_constants_plist_, "gravity");
  double gravity = gravity_plist_->get<Teuchos::Array<double>>("value")[2];

  auto perm_plist_ = Teuchos::sublist(
		  Teuchos::sublist(evals_plist_, "permeability"), "function");
  c=0;
  for (auto& entry : *perm_plist_) {
	std::string layer_name = entry.first;

    Teuchos::RCP<Teuchos::ParameterList> layer2_plist_ = Teuchos::sublist(
		                  Teuchos::sublist(
		                    Teuchos::sublist(perm_plist_, layer_name),"function"),
					    "function-constant");
    double perm = layer2_plist_->get<double>("value");
    perm = (elmdata.hksat_[c]/1000.0)*visco/den_mass/(-gravity);   //mmH2O/s --> m2
    layer2_plist_->set("value", perm);
    c++;
  }


  // (1c) water retention curve models, in 'state ->evaluators -> saturation_liquid'
  Teuchos::RCP<Teuchos::ParameterList> wrm_plist_ =
		  Teuchos::sublist(evals_plist_, "saturation_liquid");

  auto wrm_constants_plist_ = Teuchos::sublist(wrm_plist_, "WRM parameters");
  c=0;
  for (auto& entry : *wrm_constants_plist_) {
    std::string layer_name = entry.first;

    auto layer3_plist_ = Teuchos::sublist(wrm_constants_plist_, layer_name);
    std::string wrm_type = layer3_plist_->get<std::string>("WRM Type");
    if(wrm_type == "van Genuchten"){
       double vG_alpha = layer3_plist_->get<double>("van Genuchten alpha [Pa^-1]");
       double vG_n = layer3_plist_->get<double>("van Genuchten n [-]");
       double vG_sr = layer3_plist_->get<double>("residual saturation [-]");
       double vG_s0 = layer3_plist_->get<double>("smoothing interval width [saturation]", 0.0);

       // (TODO) data passing from elm

    } else if(wrm_type == "Clapp Hornberger"){
       double smpsat = layer3_plist_->get<double>("Clapp Hornberger smpsat [Pa]");
       double bsw = layer3_plist_->get<double>("Clapp Hornberger bsw [-]");
       double sr = layer3_plist_->get<double>("residual saturation [-]", 0.0);
       double s0 = layer3_plist_->get<double>("near-saturation inflection point interval [saturation]", 0.08);
       //double pcx = layer3_plist_->get<double>("dry-end smoothing starting point [Pa]", 1.0e10);

       smpsat = elmdata.CH_smpsat_[c];
       bsw = elmdata.CH_bsw_[c];
       sr = elmdata.CH_sr_[c];

       layer3_plist_->set("Clapp Hornberger smpsat [Pa]", smpsat);
       layer3_plist_->set("Clapp Hornberger bsw [-]", bsw);
       layer3_plist_->set("residual saturation [-]", sr);



    }
    c++;

  }
}

void elm_ats_plist::plist_cycle_driver_reset(const double t0, const double dt) {
  Amanzi::Utils::Units units;

  // t0_
  elm_drv_plist_->set("start time",t0,"s");
  elm_drv_plist_->set("start time units","s");

  // t1_
  elm_drv_plist_->set("end time",t0+dt,"s");
  elm_drv_plist_->set("end time units","s");

  // dt_
  elm_drv_plist_->set("max time step size [s]", dt);
  elm_drv_plist_->set("min time step size [s]", 1.0e-3);

}


} // close namespace ATS

