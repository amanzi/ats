/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Svetlana Tokareva
------------------------------------------------------------------------- */
#include "boost/algorithm/string/predicate.hpp"

#include "lake_thermo_bc_factory.hh"
#include "advection_factory.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_AdvectionUpwind.hh"
//#include "LinearOperatorFactory.hh"
#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"
#include "lake_enthalpy_evaluator.hh"
#include "density_evaluator.hh"
#include "lake_energy_evaluator.hh"
#include "thermal_conductivity_evaluator.hh"
#include "lake_heat_capacity_evaluator.hh"
#include "heat_flux_bc_evaluator.hh"
#include "lake_evaporation_rate_evaluator.hh"
#include "lake_surface_temperature_evaluator.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "lake_thermo_pk.hh"

#define MORE_DEBUG_FLAG 0


namespace Amanzi {
namespace LakeThermo {


Lake_Thermo_PK::Lake_Thermo_PK(Teuchos::ParameterList& FElist,
    const Teuchos::RCP<Teuchos::ParameterList>& plist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& solution) :
        PK(FElist, plist, S, solution),
        PK_PhysicalBDF_Default(FElist, plist, S, solution),
        modify_predictor_with_consistent_faces_(false),
        modify_predictor_for_freezing_(false),
        coupled_to_subsurface_via_temp_(false),
        coupled_to_subsurface_via_flux_(false),
        coupled_to_surface_via_temp_(false),
        coupled_to_surface_via_flux_(false),
        coupled_to_snow_(false),
        niter_(0),
        flux_exists_(true),
        implicit_advection_(true) {

  //  if (!plist_->isParameter("conserved quantity key suffix"))
  //    plist_->set("conserved quantity key suffix", "energy");
  //
  // set a default error tolerance
  if (domain_.find("surface") != std::string::npos) {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance",
        .01 * 55000.);
    soil_atol_ = 0.;
  } else {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance",
        .5 * .1 * 55000.);
    soil_atol_ = 0.5 * 2000. * 620.e-6;  // porosity * particle density soil * heat capacity soil * 1 degree
    // or, dry bulk density soil * heat capacity soil * 1 degree, in MJ
  }

  if (!plist_->isParameter("absolute error tolerance")) {
    plist_->set("absolute error tolerance", 76.e-6); // energy of 1 degree C of water per mass_atol, in MJ/mol water
  }
}


// -------------------------------------------------------------
// Setup
// -------------------------------------------------------------
void Lake_Thermo_PK::Setup() {
  PK_PhysicalBDF_Default::Setup();

  SetupLakeThermo_();
  SetupPhysicalEvaluators_();
};


void Lake_Thermo_PK::SetupLakeThermo_() {
  // Set up keys if they were not already set.
  temperature_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
  density_key_ = Keys::readKey(*plist_, domain_, "density", "density");
  energy_key_ = Keys::readKey(*plist_, domain_, "energy", "energy");
  //  wc_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  enthalpy_key_ = Keys::readKey(*plist_, domain_, "enthalpy", "enthalpy");
  flux_key_ = Keys::readKey(*plist_, domain_, "mass flux", "mass_flux");
  energy_flux_key_ = Keys::readKey(*plist_, domain_, "diffusive energy flux", "diffusive_energy_flux");
  adv_energy_flux_key_ = Keys::readKey(*plist_, domain_, "advected energy flux", "advected_energy_flux");
  conductivity_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");
  heat_capacity_key_ = Keys::readKey(*plist_, domain_, "heat capacity", "heat_capacity");
  uw_conductivity_key_ = Keys::readKey(*plist_, domain_, "upwinded thermal conductivity", "upwind_thermal_conductivity");
  cell_is_ice_key_ = Keys::readKey(*plist_, domain_, "ice", "ice");
  depth_key_ = Keys::readKey(*plist_, domain_, "water depth", "water_depth");
  surface_flux_key_ = Keys::readKey(*plist_, domain_, "surface flux", "surface_flux");
  evaporation_rate_key_ = Keys::readKey(*plist_, domain_, "evaporation rate", "evaporation_rate");
  surface_temperature_key_ = Keys::readKey(*plist_, domain_, "surface-temperature", "surface-temperature");
  std::cout << "surface_temperature_key_ = " << surface_temperature_key_ << std::endl;
  ss_temp_key_ = Keys::readKey(*plist_, domain_, "surface-temperature", "surface-temperature");

  // Check if we run a coupled simulation
  coupled_to_snow_ =
      plist_->get<bool>("coupled to snow", false);    

  std::cout << "coupled_to_snow_ = " << coupled_to_snow_ << std::endl;

  //  std::string domain_surf;
  //  domain_surf = plist_->get<std::string>("surface domain name", "surface");
  //  surface_flux_key_ = Keys::readKey(*plist_, domain_surf, "surface flux", "surface_flux");
  //  S_->RequireField(surface_flux_key_)
  //	  ->SetMesh(S_->GetMesh(domain_surf))
  //	  ->AddComponent("cell", AmanziMesh::CELL, 1);

  // Get data for special-case entities.
  S_->Require<CompositeVector, CompositeVectorSpace>(cell_vol_key_,tag_current_,name_).SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(cell_vol_key_,tag_current_);
  S_->Require<double>("atmospheric_pressure",Tags::DEFAULT);

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  LakeThermoBCFactory bc_factory(mesh_, bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_diff_flux_ = bc_factory.CreateDiffusiveFlux();
  //  bc_flux_ = bc_factory.CreateTotalFlux();

  bc_adv_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

 // -- nonlinear coefficient
  std::string method_name =
    plist_->get<std::string>("upwind conductivity method", "arithmetic mean");
  if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_, tag_next_));
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_, tag_next_));
  } else {
    Errors::Message message;
    message << "Energy PK has no upwinding method named: " << method_name;
    Exceptions::amanzi_throw(message);
  }

  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {
    S_->Require<CompositeVector, CompositeVectorSpace>(uw_conductivity_key_, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S_->Require<CompositeVector, CompositeVectorSpace>(uw_conductivity_key_, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in energy.");
    Exceptions::amanzi_throw(message);
  }
  S_->GetRecordW(uw_conductivity_key_, tag_next_, name_).set_io_vis(false);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the forward operator for the advection term
  Teuchos::ParameterList advect_plist = plist_->sublist("advection");
  matrix_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist, mesh_));
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);

  // -- create the operators for the preconditioner
  //    diffusion
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary", mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") && mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary", mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema", mfd_plist.get<Teuchos::Array<std::string> >("schema"));

  //    this needs to get some better logic now that symbolic doesn't use this. --etc
  precon_used_ = plist_->isSublist("preconditioner") ||
      plist_->isSublist("inverse") ||
      plist_->isSublist("linear solver");

  if (precon_used_) {
    auto& inv_list = mfd_pc_plist.sublist("inverse");
    inv_list.setParameters(plist_->sublist("inverse"));
    // old style... deprecate me!
    inv_list.setParameters(plist_->sublist("preconditioner"));
    inv_list.setParameters(plist_->sublist("linear solver"));
  }

  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  //    If using approximate Jacobian for the preconditioner, we also
  //    need derivative information.  This means upwinding the
  //    derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    // if (preconditioner_->RangeMap().HasComponent("face")) {
    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default"){
      // MFD -- upwind required
      dconductivity_key_ = Keys::getDerivKey(conductivity_key_, key_);
      duw_conductivity_key_ = Keys::getDerivKey(uw_conductivity_key_, key_);

      S_->Require<CompositeVector, CompositeVectorSpace>(duw_conductivity_key_, tag_next_, name_)
            .SetMesh(mesh_)->SetGhosted()
            ->SetComponent("face", AmanziMesh::FACE, 1);

      // upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
      //                                 dconductivity_key_, duw_conductivity_key_));
      upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, energy_flux_key_, 1.e-8));

    } else {
      // FV -- no upwinding
      dconductivity_key_ = Keys::getDerivKey(conductivity_key_, key_);
      duw_conductivity_key_ = "";
    }
  } else {
    dconductivity_key_ = "";
    duw_conductivity_key_ = "";
  }


  // -- accumulation terms
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  //  -- advection terms
  implicit_advection_ = !plist_->get<bool>("explicit advection", false);
  if (implicit_advection_) {
    implicit_advection_in_pc_ = !plist_->get<bool>("supress advective terms in preconditioner", false);

    if (implicit_advection_in_pc_) {
      Teuchos::ParameterList advect_plist_pc = plist_->sublist("advection preconditioner");
      preconditioner_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist_pc, preconditioner_));
      preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
    }
  }

  //  //    symbolic assemble
  //  precon_used_ = plist_->isSublist("preconditioner");
  //  if (precon_used_) {
  //    preconditioner_->SymbolicAssembleMatrix();
  //    preconditioner_->InitializePreconditioner(plist_->sublist("preconditioner"));
  //
  //    //    Potentially create a linear solver
  //    if (plist_->isSublist("linear solver")) {
  //      Teuchos::ParameterList linsolve_sublist = plist_->sublist("linear solver");
  //      AmanziSolvers::LinearOperatorFactory<Operators::Operator,CompositeVector,CompositeVectorSpace> fac;
  //      lin_solver_ = fac.Create(linsolve_sublist, preconditioner_);
  //    } else {
  //      lin_solver_ = preconditioner_;
  //    }
  //  }

  // -- advection of enthalpy
  S_->Require<CompositeVector, CompositeVectorSpace>(enthalpy_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  Teuchos::ParameterList enth_plist = plist_->sublist("enthalpy evaluator");
  std::cout << "plist_ = " << plist_ << std::endl;
  std::cout << "enth_plist = " << enth_plist << std::endl;
  // enth_plist.setParameters(plist_->sublist("enthalpy evaluator"));
  // enth_plist.set("enthalpy key", enthalpy_key_);
  enth_plist.setParameters(plist_->sublist("enthalpy evaluator"));
  enth_plist.set<std::string>("evaluator type", "enthalpy");
  std::cout << "enth_plist = " << enth_plist << std::endl;  
  Teuchos::RCP<LakeEnthalpyEvaluator> enth =
      Teuchos::rcp(new LakeEnthalpyEvaluator(enth_plist));
  std::cout << "check 1" << std::endl;
  S_->SetEvaluator(enthalpy_key_, tag_next_, enth);
  std::cout << "check 2" << std::endl;

  // -- density evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(density_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  // Teuchos::ParameterList den_plist =
  //     plist_->sublist("density evaluator");
  // den_plist.set("evaluator name", density_key_);
  // Teuchos::RCP<LakeThermo::DensityEvaluator> den =
  //     Teuchos::rcp(new LakeThermo::DensityEvaluator(den_plist));
  // S_->SetEvaluator(density_key_, den);

  // -- energy evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(energy_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  Teuchos::ParameterList enrg_plist =
      plist_->sublist("energy evaluator");
  enrg_plist.set("evaluator name", energy_key_);
  Teuchos::RCP<LakeThermo::LakeEnergyEvaluator> enrg =
      Teuchos::rcp(new LakeThermo::LakeEnergyEvaluator(enrg_plist));
  S_->SetEvaluator(energy_key_, tag_next_, enrg);

  // -- surface temperature evaluator
  std::string domain_surf_temp;
  if (domain_ == "domain" || domain_ == "") {
    domain_surf_temp = plist_->get<std::string>("surface domain name", "surface");
  } else {
    domain_surf_temp = plist_->get<std::string>("surface domain name", "surface_"+domain_);
  }
  std::cout << "domain_surf_temp = " << domain_surf_temp << std::endl;
  // S_->RequireField(surface_temperature_key_)->SetMesh(S_->GetMesh(domain_surf_temp))
  //       ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  // Teuchos::ParameterList surftemp_plist =
  //     plist_->sublist("lake surface temperature evaluator");
  // surftemp_plist.set("evaluator name", surface_temperature_key_);
  // Teuchos::RCP<LakeThermo::LakeSurfaceTemperatureEvaluator> surftemp =
  //     Teuchos::rcp(new LakeThermo::LakeSurfaceTemperatureEvaluator(surftemp_plist));
  // S_->SetEvaluator(surface_temperature_key_, surftemp);
  // std::cout << "surface_temperature_key_ = " << surface_temperature_key_ << std::endl;
  // std::cout << "after creating surface temperature evaluator" << std::endl;

  // surface temperature as independent variable
  S_->Require<CompositeVector, CompositeVectorSpace>(surface_temperature_key_,tag_next_,name_).SetMesh(S_->GetMesh(domain_surf_temp))
          ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList surftemp_plist;
  surftemp_plist.set<std::string>("evaluator name", surface_temperature_key_);
  auto surftemp = Teuchos::rcp(new EvaluatorPrimaryCV(surftemp_plist));
  AMANZI_ASSERT(S_ != Teuchos::null);
  S_->SetEvaluator(surface_temperature_key_, tag_next_, surftemp);

  S_->RequireEvaluator(surface_temperature_key_,tag_next_);

  S_->Require<CompositeVector, CompositeVectorSpace>(surface_temperature_key_, tag_next_, name_).SetMesh(S_->GetMesh(domain_surf_temp))->SetGhosted()
		      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- thermal conductivity evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(conductivity_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
      plist_->sublist("thermal conductivity evaluator");
  tcm_plist.set("evaluator name", conductivity_key_);
  Teuchos::RCP<LakeThermo::ThermalConductivityEvaluator> tcm =
      Teuchos::rcp(new LakeThermo::ThermalConductivityEvaluator(tcm_plist));
  S_->SetEvaluator(conductivity_key_, tag_next_, tcm);

  // -- heat capacity evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(heat_capacity_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList hc_plist =
      plist_->sublist("heat capacity evaluator");
  hc_plist.set("evaluator name", heat_capacity_key_);
  Teuchos::RCP<LakeThermo::LakeHeatCapacityEvaluator> hc =
      Teuchos::rcp(new LakeThermo::LakeHeatCapacityEvaluator(hc_plist));
  S_->SetEvaluator(heat_capacity_key_, tag_next_, hc);

  // -- evaporation rate evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(evaporation_rate_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList er_plist =
      plist_->sublist("evaporation rate evaluator");
  er_plist.set("evaluator name", evaporation_rate_key_);
  Teuchos::RCP<LakeThermo::LakeEvaporationRateEvaluator> er =
      Teuchos::rcp(new LakeThermo::LakeEvaporationRateEvaluator(er_plist));
  S_->SetEvaluator(evaporation_rate_key_, tag_next_, er);


  // THIS SHOULD NOT BE DEFINED ON THE MESH -- > ONLY ONE VALUE FOR THE UPPER BOUNDARY
  // -- surface flux evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(surface_flux_key_,tag_next_).SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList sflx_plist =
      plist_->sublist("surface flux evaluator");
  sflx_plist.set("evaluator name", surface_flux_key_);
  Teuchos::RCP<LakeThermo::HeatFluxBCEvaluator> sflx =
      Teuchos::rcp(new LakeThermo::HeatFluxBCEvaluator(sflx_plist));
  S_->SetEvaluator(surface_flux_key_, tag_next_, sflx);

  //  // source terms
  //  is_source_term_ = plist_->get<bool>("source term");
  //  is_source_term_differentiable_ = plist_->get<bool>("source term is differentiable", true);
  //  is_source_term_finite_differentiable_ = plist_->get<bool>("source term finite difference", false);
  //  if (is_source_term_) {
  //    source_key_ = Keys::readKey(*plist_, domain_, "source", "total_energy_source");
  //    S_->RequireField(source_key_)->SetMesh(mesh_)
  //        ->AddComponent("cell", AmanziMesh::CELL, 1);
  //    S_->RequireEvaluator(source_key_);
  //  }

  // // coupling terms
  // // -- subsurface PK, coupled to the surface
  // coupled_to_surface_via_flux_ =
  //     plist_->get<bool>("coupled to surface via flux", false);
  // if (coupled_to_surface_via_flux_) {
  //   std::string domain_surf;
  //   if (domain_ == "domain" || domain_ == "") {
  //     domain_surf = plist_->get<std::string>("surface domain name", "surface");
  //   } else {
  //     domain_surf = plist_->get<std::string>("surface domain name", "surface_"+domain_);
  //   }
  //   ss_flux_key_ = Keys::readKey(*plist_, domain_surf, "surface-subsurface energy flux", "surface_subsurface_energy_flux");
  //   S_->RequireField(ss_flux_key_)
  //           ->SetMesh(S_->GetMesh(domain_surf))
  //           ->AddComponent("cell", AmanziMesh::CELL, 1);
  // }

  // coupled_to_surface_via_temp_ =
  //     plist_->get<bool>("coupled to surface via temperature", false);
  // if (coupled_to_surface_via_temp_) {
  //   // surface temperature used for BCs
  //   S_->RequireField("surface_temperature")
  //           ->SetMesh(S_->GetMesh("surface"))
  //           ->AddComponent("cell", AmanziMesh::CELL, 1);
  // }

  // // -- Make sure coupling isn't flagged multiple ways.
  // if (coupled_to_surface_via_flux_ && coupled_to_surface_via_temp_) {
  //   Errors::Message message("Lake Thermo PK requested both flux and temperature coupling -- choose one.");
  //   Exceptions::amanzi_throw(message);
  // }

  CompositeVectorSpace matrix_cvs = matrix_->RangeMap();
  matrix_cvs.AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1); 

  // -- primary variable
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_).Update(matrix_cvs)->SetGhosted();

  // Require a field for the mass flux for advection.
  flux_exists_ = S_->HasRecord(flux_key_,tag_next_); // this bool is needed to know if PK
  // makes flux or we need an
  // independent variable evaluator

  S_->Require<CompositeVector, CompositeVectorSpace>(flux_key_, tag_next_, name_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("face", AmanziMesh::FACE, 1);

  // Require a field for water content
  //  S_->RequireField(wc_key_, name_)->SetMesh(mesh_)->SetGhosted()
  //        ->AddComponent("cell", AmanziMesh::CELL, 1);

  // Require a field for the energy fluxes.
  S_->Require<CompositeVector, CompositeVectorSpace>(energy_flux_key_, tag_next_, name_).SetMesh(mesh_)->SetGhosted()
          ->SetComponent("face", AmanziMesh::FACE, 1);
  S_->Require<CompositeVector, CompositeVectorSpace>(adv_energy_flux_key_, tag_next_, name_).SetMesh(mesh_)->SetGhosted()
          ->SetComponent("face", AmanziMesh::FACE, 1);

  // ice markers
  S_->Require<CompositeVector, CompositeVectorSpace>(cell_is_ice_key_,tag_next_,name_).SetMesh(mesh_)
          ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);

  // depth
  S_->Require<CompositeVector, CompositeVectorSpace>(depth_key_,tag_next_,name_).SetMesh(mesh_)
          ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList elist_depth;
  elist_depth.set<std::string>("evaluator name", depth_key_);
  auto eval_depth = Teuchos::rcp(new EvaluatorPrimaryCV(elist_depth));
  AMANZI_ASSERT(S_ != Teuchos::null);
  S_->SetEvaluator(depth_key_, tag_next_, eval_depth);

  S_->RequireEvaluator(depth_key_,tag_next_);

  S_->Require<CompositeVector, CompositeVectorSpace>(depth_key_, tag_next_, name_).SetMesh(mesh_)->SetGhosted()
		      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- simply limit to close to 0
  modify_predictor_for_freezing_ =
      plist_->get<bool>("modify predictor for freezing", false);
  // -- ewc and other predictors can result in odd face values
  modify_predictor_with_consistent_faces_ =
      plist_->get<bool>("modify predictor with consistent faces", false);

  // correction controls
  T_limit_ = plist_->get<double>("limit correction to temperature change [K]", -1.);
};

// -------------------------------------------------------------
// Create the physical evaluators that are specific to Lake.
// -------------------------------------------------------------
void Lake_Thermo_PK::SetupPhysicalEvaluators_() {
  // -- Density
  S_->Require<CompositeVector, CompositeVectorSpace>(density_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(density_key_,tag_next_);

  // -- Energy
  S_->Require<CompositeVector, CompositeVectorSpace>(energy_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(energy_key_,tag_next_);

  if (coupled_to_snow_) {
    // setup some evaluators needed for snow

    // -- Porosity
    // S_->RequireField("porosity")->SetMesh(mesh_)->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("porosity",tag_next_);

    // -- saturation_gas
    // S_->RequireField("saturation_gas")->SetMesh(mesh_)->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("saturation_gas",tag_next_);

    // -- saturation_gas
    // S_->RequireField("snow-precipitation")->SetMesh(S_->GetMesh("snow"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("snow-precipitation",tag_next_);

    // -- surface-air_temperature
    // S_->RequireField("surface-air_temperature")->SetMesh(S_->GetMesh("surface"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("surface-air_temperature",tag_next_);

    // -- surface-albedos
    // S_->RequireField("surface-albedos")->SetMesh(S_->GetMesh("surface"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 2);
    S_->RequireEvaluator("surface-albedos",tag_next_);

    // -- surface-ponded_depth
    // S_->RequireField("surface-ponded_depth")->SetMesh(S_->GetMesh("surface"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("surface-ponded_depth",tag_next_);

    // -- surface-mass_density_liquid
    // S_->RequireField("surface-mass_density_liquid")->SetMesh(S_->GetMesh("surface"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("surface-mass_density_liquid",tag_next_);

    // -- surface-molar_density_liquid
    // S_->RequireField("surface-molar_density_liquid")->SetMesh(S_->GetMesh("surface"))->SetGhosted()
    //         ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator("surface-molar_density_liquid",tag_next_);

  }

  // // -- Surface temperature
  // std::string domain_surf_temp;
  // if (domain_ == "domain" || domain_ == "") {
  //   domain_surf_temp = plist_->get<std::string>("surface domain name", "surface");
  // } else {
  //   domain_surf_temp = plist_->get<std::string>("surface domain name", "surface_"+domain_);
  // }
  // std::cout << "domain_surf_temp = " << domain_surf_temp << std::endl;
  // S_->RequireField(surface_temperature_key_)->SetMesh(S_->GetMesh(domain_surf_temp))->SetGhosted()
  //         ->AddComponent("cell", AmanziMesh::CELL, 1);
  // S_->RequireEvaluator(surface_temperature_key_);

  // -- Conductivity
  S_->Require<CompositeVector, CompositeVectorSpace>(conductivity_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(conductivity_key_,tag_next_);

  // -- Heat capacity
  S_->Require<CompositeVector, CompositeVectorSpace>(heat_capacity_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(heat_capacity_key_,tag_next_);

  // -- Evaporation rate
  S_->Require<CompositeVector, CompositeVectorSpace>(evaporation_rate_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(evaporation_rate_key_,tag_next_);

  // -- Enthalpy
  S_->Require<CompositeVector, CompositeVectorSpace>(enthalpy_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(enthalpy_key_,tag_next_);

  // -- Surface heat flux
  S_->Require<CompositeVector, CompositeVectorSpace>(surface_flux_key_,tag_next_).SetMesh(mesh_)->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(surface_flux_key_,tag_next_);

  std::cout << "Lake SetupPhysicalEvaluators_ DONE" << std::endl;

}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Lake_Thermo_PK::Initialize() {

  std::cout << "Lake Initialize" << std::endl;

  // initialize BDF stuff and physical domain stuff
  PK_PhysicalBDF_Default::Initialize();

  R_s_ = 0.;
  R_b_ = 0.;
  alpha_e_w_ = 0.;
  alpha_e_i_ = 0.;
  S0_ = 0.;

  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  double temp = ic_list.get<double>("initial temperature [K]");

  S_->GetW<CompositeVector>(temperature_key_, tag_next_, name_).PutScalar(temp);

  // initial depth
  h_ = ic_list.get<double>("initial depth [m]",1.5);

  // initial ice thickness
  h_ice_ = 0.;

  S_->GetW<CompositeVector>(cell_is_ice_key_, tag_next_, name_).PutScalar(false);
  S_->GetRecordW(cell_is_ice_key_, tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>(depth_key_, tag_next_, name_).PutScalar(h_);
  S_->GetRecordW(depth_key_, tag_next_, name_).set_initialized();

#if MORE_DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << i;
    S_->GetW<CompositeVector>(namestream.str(),tag_next_,name_).PutScalar(0.);
    S_->GetRecordW(namestream.str(),tag_next_,name_).set_initialized();

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << i;
    S_->GetW<CompositeVector>(solnstream.str(),tag_next_,name_).PutScalar(0.);
    S_->GetRecordW(solnstream.str(),tag_next_,name_).set_initialized();
  }

#endif

  // initialize energy flux
  S_->GetW<CompositeVector>(energy_flux_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(energy_flux_key_, tag_next_, name_).set_initialized();
  S_->GetW<CompositeVector>(adv_energy_flux_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(adv_energy_flux_key_, tag_next_, name_).set_initialized();
  S_->GetW<CompositeVector>(uw_conductivity_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(uw_conductivity_key_, tag_next_, name_).set_initialized();
  if (!duw_conductivity_key_.empty()) {
    S_->GetW<CompositeVector>(duw_conductivity_key_, tag_next_, name_).PutScalar(0.);
    S_->GetRecordW(duw_conductivity_key_, tag_next_, name_).set_initialized();
  }

  UpdateConductivityData_(tag_current_);
  Teuchos::RCP<const CompositeVector> conductivity =
      S_->GetPtrW<CompositeVector>(uw_conductivity_key_,tag_next_,name_);

  const Epetra_MultiVector& uw_cond_c = *S_->Get<CompositeVector>(uw_conductivity_key_,tag_next_).ViewComponent("face", false);
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);

  // potentially initialize mass flux
  //  if (!flux_exists_) {
  //    S_->GetField(flux_key_, name_)->Initialize(plist_->sublist(flux_key_));
  //  }
  S_->GetW<CompositeVector>(flux_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(flux_key_, tag_next_, name_).set_initialized();

  //  S_->GetW<CompositeVector>(wc_key_, name_).PutScalar(1.0);
  //  S_->GetRecordW(wc_key_, name_).set_initialized();

  /*
  // get temperature
  const Epetra_MultiVector& temp = *S_->GetW<CompositeVector>(temperature_key_, tag_next_)
        .ViewComponent("cell",false);

  double T0 = 280., T1 = 400.;
  double zm = 0.5, Tm = 290.;
  double d = T0;
  double b = (Tm-T0-zm*zm*(T1-T0))/(zm*(1.-zm));
  double a = T1-T0-b;

  int ncells_owned = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      temp[0][c] = a*xc[2]*xc[2] + b*xc[2] + d;
  }
   */

  S_->GetRecordW(temperature_key_, tag_next_, name_).set_initialized();

  int ncomp = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  const Epetra_MultiVector& temp_v = *S_->Get<CompositeVector>(temperature_key_,tag_next_).ViewComponent("cell",false);
  std::cout << "surface_temperature_key_ = " << surface_temperature_key_ << std::endl;
  Epetra_MultiVector surf_temp = *S_->Get<CompositeVector>(surface_temperature_key_,tag_next_).ViewComponent("cell",false);
  std::cout << "name_ = " << name_ << std::endl;
  S_->GetW<CompositeVector>(surface_temperature_key_, tag_next_, name_).PutScalar(temp_v[0][ncomp-1]);

  S_->GetRecordW(surface_temperature_key_, tag_next_, name_).set_initialized();

  // summary of initialization
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Lake Thermo PK initialized." << std::endl;

  std::cout << "Lake Initialize DONE" << std::endl;  

}


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  double dt = t_new - t_old;
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, tag_next);

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S_->GetPtr<CompositeVector>(temperature_key_,tag_next);
  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);

  // cell volumes
  const Epetra_MultiVector& cv =
         *S_->Get<CompositeVector>(Keys::getKey(domain_,"cell_volume"),tag_next).ViewComponent("cell",false);

  // ice marker
  const Epetra_MultiVector& ice =
         *S_->Get<CompositeVector>(cell_is_ice_key_,tag_next).ViewComponent("cell",false);
  S_->GetW<CompositeVector>(cell_is_ice_key_, tag_next_, name_).PutScalar(false);

  int ncomp = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);

  Epetra_MultiVector surf_temp = *S_->Get<CompositeVector>(surface_temperature_key_,tag_next).ViewComponent("cell",false);
  std::cout << "name_ = " << name_ << std::endl;
  S_->GetW<CompositeVector>(surface_temperature_key_, tag_next, name_).PutScalar(temp_v[0][ncomp-1]);
  surf_temp[0][0] = temp_v[0][ncomp-1];

  std::cout << "surface_temperature_key_ = " << surface_temperature_key_ << std::endl;
  std::cout << "Updated surface temperature to " << surf_temp[0][0] << std::endl;


  bool ice_cover_ = false; // first always assume that there is no ice

  i_ice_max = 0;
  i_ice_min = ncomp-1;
  d_ice = 0;
  h_ice_ = 0;

  for (int i=0; i!=ncomp; ++i) {
    if (temp_v[0][i] < 273.15) { // check if there is ice cover
      ice_cover_ = true;
      i_ice_max = i;
      ice[0][i] = true;
      d_ice++;
    }
  } // i

  for (int i=ncomp-1; i!=-1; --i) {
    if (temp_v[0][i] < 273.15) { // check if there is ice cover
      ice_cover_ = true;
      i_ice_min = i;
    }
  } // i

  int d_thawed = ncomp-1-i_ice_max;  // thickness of thawed layer [cells]
//  int d_ice = i_ice_max-i_ice_min+1; // thickness of ice layer [cells]

  h_ice_prev = h_ice_;
  // if (d_ice > 0) h_ice_ = d_ice*cv[0][0]*h_; // in [m], assume uniform mesh

  std::vector<double> temp_new(ncomp); // new temperatures for swapping the cells

  for (int i=0; i!=ncomp; ++i) {
    temp_new[i] = temp_v[0][i]; //-100.; //temp_v[0][i];
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
      std::cout << "copy cell " << i_ice_max+1+i << " to " << i_ice_min+i << std::endl;
    } // i

    for (int i=0; i!=ncomp; ++i) {
    //  temp_v[0][i] = temp_new[i];
    }

    std::cout << "Temperature after swap " << std::endl;
    for (int i=ncomp-1; i!=-1; --i) {
      std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
      if (temp_v[0][i] < 0.) exit(0);
    }

//    exit(0);
  }

  i_ice_max = 0;
  i_ice_min = ncomp-1;
  d_ice = 0;
  h_ice_ = 0;

  for (int i=0; i!=ncomp; ++i) {
    if (temp_v[0][i] < 273.15) { // check if there is ice cover
      ice_cover_ = true;
      i_ice_max = i;
      ice[0][i] = true;
      d_ice++;
    }
  } // i

  for (int i=ncomp-1; i!=-1; --i) {
    if (temp_v[0][i] < 273.15) { // check if there is ice cover
      ice_cover_ = true;
      i_ice_min = i;
    }
  } // i

  d_thawed = ncomp-1-i_ice_max;  // thickness of thawed layer [cells]
//  d_ice = i_ice_max-i_ice_min+1; // thickness of ice layer [cells]

  h_ice_prev = h_ice_;
  // if (d_ice > 0) h_ice_ = d_ice*cv[0][0]*h_; // in [m], assume uniform mesh

  Teuchos::ParameterList& param_list = plist_->sublist("met data");
  FunctionFactory fac;
  Teuchos::RCP<Function> r_func_ = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
  Teuchos::RCP<Function> E_func_ = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));
  Teuchos::RCP<Function> SS_func_ = Teuchos::rcp(fac.Create(param_list.sublist("solar radiation")));
  Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("atmospheric downward radiation")));
  Teuchos::RCP<Amanzi::Function> T_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("air temperature")));

  std::vector<double> args(1);
  args[0] = S_->get_time();
  r_ = (*r_func_)(args);
  E_ = (*E_func_)(args);
  double SS = (*SS_func_)(args);
  double E_a = (*E_a_func_)(args);
  double T_a = (*T_a_func_)(args);

  // Compute evaporartion rate
  S_->GetEvaluator(evaporation_rate_key_,tag_next).Update(*S_, name_);
  const Epetra_MultiVector& E_v =
      *S_->Get<CompositeVector>(evaporation_rate_key_,tag_next).ViewComponent("cell",false);
  E_ = E_v[0][0]; // same everywhere

  // // set precipitation and evaporation to zero in the winter (for now because we don't have a snow layer anyway)
  // r_ = (temp_v[0][ncomp-1] < 273.15) ? 0. : r_;
  // E_ = (temp_v[0][ncomp-1] < 273.15) ? 0. : E_;

  // update depth
  double dhdt = r_ - E_ - R_s_ - R_b_;

  dhdt = (temp_v[0][ncomp-1] < 273.15) ? 0. : dhdt;

  h_ += dhdt*dt;

  if (d_ice > 0) h_ice_ = double(d_ice)*cv[0][0]*h_; ///1.5; //*1.5; //h_; //h_/ncomp;  //*cv[0][0]/h_; // in [m], assume uniform mesh

  double tpl_kappa_w  = 5.46e-1; // Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  double tpl_rho_I    = 9.1e2;   // Density of ice [kg m^{-3}]
  double tpl_L_f      = 3.3e5;   // Latent heat of fusion [J kg^{-1}]
  double Phi_T_pr0_1  = 40./3.;  // Constant in the expression for the T shape-function derivative 
  double Phi_T_pr0_2  = 20./3.;  // Constant in the expression for the T shape-function derivative 
  double depth_w = h_;
  double T_bot_p_flk = temp_v[0][0];
  double T_wML_p_flk = temp_v[0][ncomp-1];
  double Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w;
  double C_T_p_flk = 0.5;
  double Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2;         // d\Phi(0)/d\zeta (thermocline)
  Q_w_flk = Q_w_flk*fmax(Phi_T_pr0_flk, 1.);           // Account for an increased d\Phi(0)/d\zeta 
  double d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f;
  // d_h_ice_dt = (temp_v[0][ncomp-1] > 273.15) ? 0. : d_h_ice_dt;
  // h_ice_ += d_h_ice_dt*dt;

  // compute freeze rate
  double freeze_rate = (h_ice_-h_ice_prev)/dt*86400./2.54*100.;
  if (freeze_rate > 0.) std::cout << "freeze rate = " << freeze_rate << " inch/day" << std::endl;

  S_->GetEvaluator(surface_flux_key_,tag_next).Update(*S_, name_);
  const Epetra_MultiVector& flux =
      *S_->Get<CompositeVector>(surface_flux_key_,tag_next).ViewComponent("cell",false);

  // get temperature
  const Epetra_MultiVector& cond_v = *S_->Get<CompositeVector>(conductivity_key_,tag_next)
      .ViewComponent("cell",false);

  double bc_flux = flux[0][ncomp-1]/cond_v[0][ncomp-1]*h_; ///cv[0][ncomp-1];

  // save file with depth
  if (int(t_new) % 86400 == 0) {
    std::string ncells(std::to_string(ncomp));

    std::ofstream temp_distr;
    temp_distr.open ("temperature_distr_"+ncells+".txt", std::ios::out);
    for (int i = 0; i < ncomp; i++) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(i);
      temp_distr << xc[2] << " " << temp_v[0][i] << "\n";
    }

    std::ofstream tempfile;
    tempfile.open ("depth_lake_"+ncells+".txt", std::ios::app);
    tempfile << int(t_new) / 86400 << " " << h_ << " ";
    tempfile << "\n";

    std::ofstream tempfile_ice;
    tempfile_ice.open ("depth_ice_"+ncells+".txt", std::ios::app);
    tempfile_ice << int(t_new) / 86400 << " " << h_ice_ << " " << d_ice << " " << i_ice_min << " " << i_ice_max;
    tempfile_ice << "\n";

    if (coupled_to_snow_) {
      const Epetra_MultiVector& snow_depth_v = *S_->Get<CompositeVector>("snow-depth",tag_next).ViewComponent("cell",false);
      std::ofstream tempfile_snow;
      tempfile_snow.open ("depth_snow_"+ncells+".txt", std::ios::app);
      tempfile_snow << int(t_new) / 86400 << " " << snow_depth_v[0][0] << " ";
      tempfile_snow << "\n";
    }

    std::ofstream tempfile_freeze;
    tempfile_freeze.open ("freeze_rate_"+ncells+".txt", std::ios::app);
    tempfile_freeze << int(t_new) / 86400 << " " << freeze_rate << " ";
    tempfile_freeze << "\n";

    std::ofstream tempfile_precip;
    tempfile_precip.open ("precip_rate_"+ncells+".txt", std::ios::app);
    tempfile_precip << int(t_new) / 86400 << " " << r_ << " ";
    tempfile_precip << "\n";

    std::ofstream tempfile_evap;
    tempfile_evap.open ("evap_rate_"+ncells+".txt", std::ios::app);
    tempfile_evap << int(t_new) / 86400 << " " << E_ << " ";
    tempfile_evap << "\n";

    std::ofstream tempfile_sol_rad;
    tempfile_sol_rad.open ("solar_rad_"+ncells+".txt", std::ios::app);
    tempfile_sol_rad << int(t_new) / 86400 << " " << SS << " ";
    tempfile_sol_rad << "\n";

    std::ofstream tempfile_atm_down_rad;
    tempfile_atm_down_rad.open ("atm_down_rad_"+ncells+".txt", std::ios::app);
    tempfile_atm_down_rad << int(t_new) / 86400 << " " << E_a << " ";
    tempfile_atm_down_rad << "\n";

    std::ofstream tempfile_surftemp;
    tempfile_surftemp.open ("surf_temp_"+ncells+".txt", std::ios::app);
    tempfile_surftemp << int(t_new) / 86400 << " " << temp_v[0][ncomp-1] << " ";
    tempfile_surftemp << "\n";

    std::ofstream tempfile_ice_water_temp;
    tempfile_ice_water_temp.open ("ice_water_interface_temp_"+ncells+".txt", std::ios::app);
    tempfile_ice_water_temp << int(t_new) / 86400 << " " << temp_v[0][std::min(ncomp-1,i_ice_min)] << " ";
    tempfile_ice_water_temp << "\n";

    std::ofstream tempfile_bottomtemp;
    tempfile_bottomtemp.open ("bottom_temp_"+ncells+".txt", std::ios::app);
    tempfile_bottomtemp << int(t_new) / 86400 << " " << temp_v[0][0] << " ";
    tempfile_bottomtemp << "\n";

    std::ofstream tempfile_surfflux;
    tempfile_surfflux.open ("surf_flux_"+ncells+".txt", std::ios::app);
    tempfile_surfflux << int(t_new) / 86400 << " " << bc_flux << " ";
    tempfile_surfflux << "\n";

    std::ofstream tempfile_airtemp;
    tempfile_airtemp.open ("air_temp_"+ncells+".txt", std::ios::app);
    tempfile_airtemp << int(t_new) / 86400 << " " << T_a << " ";
    tempfile_airtemp << "\n";

  }

  // create a mesh
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, plist_->sublist("regions"), *comm));
  bool request_faces = true, request_edges = true;
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Amanzi::AmanziMesh::Preference({Amanzi::AmanziMesh::Framework::MSTK}));

  std::cout << "h_ = " << h_ << std::endl;
//  mesh_scaled_ = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, h_, 1, 1, ncomp, request_faces, request_edges);

  mesh_ = S_->GetMesh(domain_);
  std::cout << "domain_ = " << domain_ << std::endl;
  std::cout << "mesh_ = " << mesh_ << std::endl;
  // setVisMesh(mesh_);
//  mesh_->vis_mesh_ = mesh_; //mesh_scaled_;

  bc_temperature_->Compute(S_->get_time());
  bc_diff_flux_->Compute(S_->get_time());
  //  bc_flux_->Compute(S_->get_time());
  UpdateBoundaryConditions_(tag_next);

  niter_ = 0;
  bool update = UpdateConductivityData_(tag_next);
  update |= S_->GetEvaluator(key_,tag_next).Update(*S_, name_);

  if (update) {
    Teuchos::RCP<const CompositeVector> temp = S_->GetPtr<CompositeVector>(key_,tag_next);
    Teuchos::RCP<const CompositeVector> conductivity = S_->GetPtr<CompositeVector>(uw_conductivity_key_,tag_next);
    matrix_diff_->global_operator()->Init();
    matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
    matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
    matrix_diff_->ApplyBCs(true, true, true);

    Teuchos::RCP<CompositeVector> eflux = S_->GetPtrW<CompositeVector>(energy_flux_key_, tag_next, name_);
    matrix_diff_->UpdateFlux(temp.ptr(), eflux.ptr());

    // calculate the advected energy as a diagnostic
    Teuchos::RCP<const CompositeVector> flux = S_->GetPtr<CompositeVector>(flux_key_,tag_next);
    matrix_adv_->Setup(*flux);
    S_->GetEvaluator(enthalpy_key_,tag_next).Update(*S_, name_);
    Teuchos::RCP<const CompositeVector> enth = S_->GetPtr<CompositeVector>(enthalpy_key_,tag_next);
    ApplyDirichletBCsToTemperature_(tag_next);

    Teuchos::RCP<CompositeVector> adv_energy = S_->GetPtrW<CompositeVector>(adv_energy_flux_key_, tag_next, name_);  
    matrix_adv_->UpdateFlux(enth.ptr(), flux.ptr(), bc_adv_, adv_energy.ptr());
  }
}

/* ******************************************************************
 * Returns the first cell attached to a boundary face.
 ****************************************************************** */
int Lake_Thermo_PK::BoundaryFaceGetCell(int f) const
{
  auto cells = mesh_->getFaceCells(f);
  return cells[0];
}

/* ******************************************************************
 * For FV methods, set the boundary face temp.
 ****************************************************************** */
void Lake_Thermo_PK::ApplyDirichletBCsToBoundaryFace_(const Tag& tag)
{

  auto& T_vec = *S_->GetPtrW<CompositeVector>(key_, tag, name_);
  Epetra_MultiVector& temp_bf = *T_vec.ViewComponent("boundary_face", false);
  const Epetra_MultiVector& temp_c = *T_vec.ViewComponent("cell", false);

  const std::vector<int>& bc_model = bc_->bc_model();
  const std::vector<double>& bc_value = bc_->bc_value();

  int nbfaces = temp_bf.MyLength();

  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      temp_bf[0][bf] = bc_value[f];
    } else {
      temp_bf[0][bf] = temp_c[0][BoundaryFaceGetCell(f)];
    }
  }
}


bool Lake_Thermo_PK::UpdateConductivityData_(const Tag& tag) {
  Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(conductivity_key_, tag);
  bool update = S_->GetEvaluator(conductivity_key_,tag).Update(*S_, name_);
  if (update) {

    Teuchos::RCP<CompositeVector> uw_cond =
      S_->GetPtrW<CompositeVector>(uw_conductivity_key_, tag, name_);
    upwinding_->Update(*cond, *uw_cond, *S_);

     if (uw_cond->HasComponent("face")) uw_cond->ScatterMasterToGhosted("face");
  }
  return update;
}


bool Lake_Thermo_PK::UpdateConductivityDerivativeData_(const Tag& tag) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating conductivity derivatives? ";

  bool update = S_->GetEvaluator(conductivity_key_,tag).UpdateDerivative(*S_, name_, key_, tag);

  if (update) {
    Teuchos::RCP<const CompositeVector> dcond =
      S_->GetDerivativePtr<CompositeVector>(conductivity_key_, tag, key_, tag);
    if (!duw_conductivity_key_.empty()) {
      // get upwind conductivity data
      Teuchos::RCP<CompositeVector> duw_cond =
        S_->GetPtrW<CompositeVector>(duw_conductivity_key_, tag, name_);
      duw_cond->PutScalar(0.);
      upwinding_deriv_->Update(*dcond, *duw_cond, *S_); 
    } else {
      Teuchos::RCP<const CompositeVector> dcond =
          S_->GetPtr<CompositeVector>(dconductivity_key_,tag);
      dcond->ScatterMasterToGhosted("cell");
    }
  }
  return update;
}

// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::UpdateBoundaryConditions_(const Tag& tag) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  auto& adv_markers = bc_adv_->bc_model();
  auto& adv_values = bc_adv_->bc_value();

  for (unsigned int n=0; n!=markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
    adv_markers[n] = Operators::OPERATOR_BC_NONE;
    adv_values[n] = 0.0;
  }

  bool snow_exists = false;

  if (coupled_to_snow_) {
    // get snow temperature
    const Epetra_MultiVector& snow_temp_v = *S_->Get<CompositeVector>("snow-temperature",tag).ViewComponent("cell",false);
    std::cout << "T_snow = " << snow_temp_v[0][0] << std::endl;
    const Epetra_MultiVector& snow_depth_v = *S_->Get<CompositeVector>("snow-depth",tag).ViewComponent("cell",false);
    std::cout << "h_snow = " << snow_depth_v[0][0] << std::endl;
    const Epetra_MultiVector& surf_temp_v = *S_->Get<CompositeVector>("surface-temperature",tag).ViewComponent("cell",false);
    std::cout << "T_surf = " << surf_temp_v[0][0] << std::endl;
    const Epetra_MultiVector& surf_qE_cond_v = *S_->Get<CompositeVector>("surface-qE_conducted",tag).ViewComponent("cell",false);
    std::cout << "surf_qE_cond_v = " << surf_qE_cond_v[0][0] << std::endl; 
    if (snow_depth_v[0][0] > 2.e-2 && h_ice_ > 0.) snow_exists = true;
  }

  if (snow_exists) { // snow exists 

    // get snow temperature
    const Epetra_MultiVector& snow_temp_v = *S_->Get<CompositeVector>("snow-temperature",tag).ViewComponent("cell",false);
    std::cout << "T_snow = " << snow_temp_v[0][0] << std::endl;
    const Epetra_MultiVector& snow_depth_v = *S_->Get<CompositeVector>("snow-depth",tag).ViewComponent("cell",false);
    std::cout << "h_snow = " << snow_depth_v[0][0] << std::endl;
    const Epetra_MultiVector& surf_temp_v = *S_->Get<CompositeVector>("surface-temperature",tag).ViewComponent("cell",false);
    std::cout << "T_surf = " << surf_temp_v[0][0] << std::endl;
    const Epetra_MultiVector& surf_qE_cond_v = *S_->Get<CompositeVector>("surface-qE_conducted",tag).ViewComponent("cell",false);
    std::cout << "surf_qE_cond_v = " << surf_qE_cond_v[0][0] << std::endl; 

    std::cout << "case 1: snow exists" << std::endl;

    const Epetra_MultiVector& cv =
        *S_->Get<CompositeVector>(Keys::getKey(domain_,"cell_volume"),tag).ViewComponent("cell",false);

    int ncomp = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);

    // Dirichlet temperature boundary conditions
    for (Functions::BoundaryFunction::Iterator bc=bc_temperature_->begin();
        bc!=bc_temperature_->end(); ++bc) {
      AmanziMesh::Entity_ID_List cells;

      int f = bc->first;
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      auto fcells = mesh_->getFaceCells(f);
      values[f] = bc->second; 
      adv_values[f] = bc->second; 
      // if (fcells[0] == 0) { //bottom
      //   values[f] = bc->second; 
      //   adv_values[f] = bc->second; 
      // } else { // top
      //   double Ks = 0.029;
      //   values[f] = surf_qE_cond_v[0][0]*snow_depth_v[0][0]/Ks + snow_temp_v[0][0];
      //   std::cout << "T_bc = " << values[f] << std::endl;
      //   adv_values[f] = surf_qE_cond_v[0][0]*snow_depth_v[0][0]/Ks + snow_temp_v[0][0]; //bc->first;
      // }
    }

    // Neumann diffusive flux, not Neumann TOTAL flux.  Potentially advective flux.
    for (Functions::BoundaryFunction::Iterator bc=bc_diff_flux_->begin();
        bc!=bc_diff_flux_->end(); ++bc) {
      int f = bc->first;
     auto fcells = mesh_->getFaceCells(f);
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      // if (fcells[0] == 0) { //bottom
      //   values[f] = 0.;
      // } else {
        values[f] = -surf_qE_cond_v[0][0]/h_; ///cv[0][ncomp-1]; 
      // }
      adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
      adv_values[f] = 0.;
    }

  } else { // no snow

    std::cout << "case 2: no snow" << std::endl;

    //  // Neumann flux boundary conditions
    //  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
    //       bc!=bc_flux_->end(); ++bc) {
    //    int f = bc->first;
    //    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    //    values[f] = bc->second;
    //    adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
    //    // push all onto diffusion, assuming that the incoming enthalpy is 0 (likely mass flux is 0)
    //  }

    // Dirichlet temperature boundary conditions
    for (Functions::BoundaryFunction::Iterator bc=bc_temperature_->begin();
        bc!=bc_temperature_->end(); ++bc) {
      AmanziMesh::Entity_ID_List cells;

      int f = bc->first;
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      auto fcells = mesh_->getFaceCells(f);
      values[f] = bc->second; 
      adv_values[f] = bc->second; 
    }

  std::cout << "in BC before flux" << std::endl;
    S_->GetEvaluator(surface_flux_key_,tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& flux =
        *S_->Get<CompositeVector>(surface_flux_key_,tag).ViewComponent("cell",false);
  std::cout << "in BC after flux" << std::endl;


    // get conductivity
    const Epetra_MultiVector& cond_v = *S_->Get<CompositeVector>(conductivity_key_,tag)
        .ViewComponent("cell",false);

    const Epetra_MultiVector& cv =
            *S_->Get<CompositeVector>(Keys::getKey(domain_,"cell_volume"),tag).ViewComponent("cell",false);

    int ncomp = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
    // Neumann diffusive flux, not Neumann TOTAL flux.  Potentially advective flux.
    for (Functions::BoundaryFunction::Iterator bc=bc_diff_flux_->begin();
        bc!=bc_diff_flux_->end(); ++bc) {
      int f = bc->first;
      auto fcells = mesh_->getFaceCells(f);
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      // if (fcells[0] == 0) { //bottom
      //   values[f] = 0.;
      // } else {
        // values[f] = flux[0][ncomp-1]*h_; 
        // values[f] = flux[0][ncomp-1]*h_/cv[0][ncomp-1]; 

        // values[f] = -surf_qE_cond_v[0][0]; ///h_;

        values[f] = flux[0][ncomp-1]/h_; ///cv[0][ncomp-1]; // /cv[0][ncomp-1]; ///cond_v[0][ncomp-1];

        // values[f] = flux[0][ncomp-1]/cond_v[0][ncomp-1]*h_;
        // values[f] = flux[0][ncomp-1]/cond_v[0][ncomp-1]*h_/cv[0][ncomp-1];
      // }
      adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
      adv_values[f] = 0.; //flux[0][ncomp-1]/h_; ///cv[0][ncomp-1]; 
    }

  }

  // // Dirichlet temperature boundary conditions from a coupled surface.
  // if (coupled_to_surface_via_temp_) {
  //   // Face is Dirichlet with value of surface temp
  //   Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh(Keys::getDomain(ss_flux_key_));
  //   const Epetra_MultiVector& temp = *S_->GetFieldData("surface_temperature")
  //           ->ViewComponent("cell",false);

  //   int ncells_surface = temp.MyLength();
  //   for (int c=0; c!=ncells_surface; ++c) {
  //     // -- get the surface cell's equivalent subsurface face
  //     AmanziMesh::Entity_ID f =
  //         surface->entity_get_parent(AmanziMesh::CELL, c);

  //     // -- set that value to dirichlet
  //     markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     values[f] = temp[0][c];
  //     adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     adv_values[f] = temp[0][c];
  //   }
  // }

  //  // surface coupling
  //  if (coupled_to_surface_via_flux_) {
  //    // Diffusive fluxes are given by the residual of the surface equation.
  //    // Advective fluxes are given by the surface temperature and whatever flux we have.
  //    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh(Keys::getDomain(ss_flux_key_));
  //    const Epetra_MultiVector& flux =
  //      *S_->GetFieldData(ss_flux_key_)->ViewComponent("cell",false);
  //
  //    int ncells_surface = flux.MyLength();
  //    for (int c=0; c!=ncells_surface; ++c) {
  //      // -- get the surface cell's equivalent subsurface face
  //      AmanziMesh::Entity_ID f =
  //        surface->entity_get_parent(AmanziMesh::CELL, c);
  //
  //      // -- set that value to Neumann
  //      markers[f] = Operators::OPERATOR_BC_NEUMANN;
  //      // flux provided by the coupler is in units of J / s, whereas Neumann BCs are J/s/A
  //      values[f] = flux[0][c] / mesh_->face_area(f);
  //
  //      // -- mark advective BCs as Dirichlet: this ensures the surface
  //      //    temperature is picked up and advective fluxes are treated
  //      //    via advection operator, not diffusion operator.
  ////      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  //    }
  //  }

  // mark all remaining boundary conditions as zero diffusive flux conditions
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      auto cells = mesh_->getFaceCells(f);
      int ncells = cells.size();

      if (ncells == 1) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
        adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
        adv_values[f] = 0.;
      }
    }
  }

  // set the face temperature on boundary faces
  auto& temp = *S_->GetPtrW<CompositeVector>(key_, tag, name_);
  if (temp.HasComponent("boundary_face")) {
    ApplyDirichletBCsToBoundaryFace_(tag);
  }
};


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool Lake_Thermo_PK::IsAdmissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->Data();
  double minT, maxT;

  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  double minT_c(1.e6), maxT_c(-1.e6);
  int min_c(-1), max_c(-1);
  for (int c=0; c!=temp_c.MyLength(); ++c) {
    if (temp_c[0][c] < minT_c) {
      minT_c = temp_c[0][c];
      min_c = c;
    }
    if (temp_c[0][c] > maxT_c) {
      maxT_c = temp_c[0][c];
      max_c = c;
    }
  }

  double minT_f(1.e6), maxT_f(-1.e6);
  int min_f(-1), max_f(-1);
  if (temp->HasComponent("face")) {
    const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
    for (int f=0; f!=temp_f.MyLength(); ++f) {
      if (temp_f[0][f] < minT_f) {
        minT_f = temp_f[0][f];
        min_f = f;
      }
      if (temp_f[0][f] > maxT_f) {
        maxT_f = temp_f[0][f];
        max_f = f;
      }
    }
    minT = std::min(minT_c, minT_f);
    maxT = std::max(maxT_c, maxT_f);

  } else {
    minT = minT_c;
    maxT = maxT_c;
  }

  double minT_l = minT;
  double maxT_l = maxT;
  mesh_->getComm()->MaxAll(&maxT_l, &maxT, 1);
  mesh_->getComm()->MinAll(&minT_l, &minT, 1);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
      Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();



  if (minT < 0.0 || maxT > 10000000.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:" << std::endl;
      ENorm_t global_minT_c, local_minT_c;
      ENorm_t global_maxT_c, local_maxT_c;

      local_minT_c.value = minT_c;
      local_minT_c.gid = temp_c.Map().GID(min_c);
      local_maxT_c.value = maxT_c;
      local_maxT_c.gid = temp_c.Map().GID(max_c);

      MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

      *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
          << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

      if (temp->HasComponent("face")) {
        const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
        ENorm_t global_minT_f, local_minT_f;
        ENorm_t global_maxT_f, local_maxT_f;

        local_minT_f.value = minT_f;
        local_minT_f.gid = temp_f.Map().GID(min_f);
        local_maxT_f.value = maxT_f;
        local_maxT_f.gid = temp_f.Map().GID(max_f);


        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        *vo_->os() << "   cells (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
            << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
      }
    }
    return false;
  }
  return true;
}


// -----------------------------------------------------------------------------
// BDF takes a prediction step -- make sure it is physical and otherwise ok.
// -----------------------------------------------------------------------------
bool Lake_Thermo_PK::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
    Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  bc_temperature_->Compute(S_->get_time(tag_next_));
  bc_diff_flux_->Compute(S_->get_time(tag_next_));
  //  bc_flux_->Compute(S_next_->get_time());
  UpdateBoundaryConditions_(tag_next_);

  //  // push Dirichlet data into predictor
  //  if (u->Data()->HasComponent("boundary_cell")) {
  //    ApplyBoundaryConditions_(u->Data().ptr());
  //  }

  bool modified = false;
  if (modify_predictor_for_freezing_) {
    const Epetra_MultiVector& u0_c = *u0->Data()->ViewComponent("cell",false);
    Epetra_MultiVector& u_c = *u->Data()->ViewComponent("cell",false);

    for (int c=0; c!=u0_c.MyLength(); ++c) {
      if (u0_c[0][c] > 273.15 && u_c[0][c] < 273.15) {
        u_c[0][c] = 273.15 - .00001;
        modified = true;
      }
    }
  }

  if (modify_predictor_with_consistent_faces_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  modifications for consistent face temperatures." << std::endl;
    CalculateConsistentFaces(u->Data().ptr());
    modified = true;
  }
  return modified;
}


// -----------------------------------------------------------------------------
// Given an arbitrary set of cell values, calculate consitent face constraints.
//
//  This is useful for prediction steps, hacky preconditioners, etc.
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {

  // average cells to faces to give a reasonable initial guess
  u->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell",true);
  Epetra_MultiVector& u_f = *u->ViewComponent("face",false);

  int f_owned = u_f.MyLength();
  for (int f=0; f!=f_owned; ++f) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += u_c[0][cells[n]];
    }
    u_f[0][f] = face_value / ncells;
  }
  ChangedSolution();

  // use old BCs
  // // update boundary conditions
  // bc_temperature_->Compute(S_next_->get_time());
  // bc_diff_flux_->Compute(S_next_->get_time());
  // bc_flux_->Compute(S_next_->get_time());
  // UpdateBoundaryConditions_(S_next_.ptr());

  // use old conductivity
  // div K_e grad u
  //  bool update = UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
       S_->GetPtr<CompositeVector>(uw_conductivity_key_,tag_next_);

  // Update the preconditioner
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, u);
  //matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  matrix_diff_->ApplyBCs(true, true, true);

  // derive the consistent faces, involves a solve
  matrix_diff_->UpdateConsistentFaces(*u);
}



AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Lake_Thermo_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> du) {

  int my_limited = 0;
  int n_limited = 0;
  if (T_limit_ > 0.) {
    for (CompositeVector::name_iterator comp=du->Data()->begin();
        comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max temperature correction (" << *comp << ") = " << max << std::endl;
      }

      for (int c=0; c!=du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > T_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * T_limit_;
          my_limited++;
        }
      }
    }
    mesh_->getComm()->SumAll(&my_limited, &n_limited, 1);
  }

  if (n_limited > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limited by temperature." << std::endl;
    }
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }
  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

/* ******************************************************************
 * Return a pointer to a local operator
 ****************************************************************** */
Teuchos::RCP<Operators::Operator> Lake_Thermo_PK::my_operator(
    const Operators::OperatorType& type)
{
  if (type == Operators::OPERATOR_MATRIX) return matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW) return preconditioner_;
  return Teuchos::null;
}

} // namespace LakeThermo
} // namespace Amanzi
