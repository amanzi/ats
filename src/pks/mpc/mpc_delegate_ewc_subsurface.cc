/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.
------------------------------------------------------------------------- */
#define DEBUG_FLAG 1

#define EWC_THAWING 1
#define EWC_SATURATION 1
#define EWC_INCREASING_PRESSURE 1

#define EWC_PC_THAWING 0
#define EWC_PC_SATURATION 0
#define EWC_PC_INCREASING_PRESSURE 0

#include "ewc_model.hh"
#include "Evaluator.hh"
#include "mpc_delegate_ewc_subsurface.hh"


namespace Amanzi {

MPCDelegateEWCSubsurface::MPCDelegateEWCSubsurface(Teuchos::ParameterList& plist,
                                                   const Teuchos::RCP<State>& S)
  : MPCDelegateEWC(plist, S)
{}


bool
MPCDelegateEWCSubsurface::modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell", false);

  const double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "  Modifying predictor using SmartEWC algorithm" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("p_extrap");
    vnames.push_back("T_extrap");
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(pres_guess.ptr());
    vecs.push_back(temp_guess.ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  // T, p at the previous step
  const Epetra_MultiVector& T1 =
    *S_->GetPtr<CompositeVector>(temp_key_, tag_current_)->ViewComponent("cell", false);
  const Epetra_MultiVector& p1 =
    *S_->GetPtr<CompositeVector>(pres_key_, tag_current_)->ViewComponent("cell", false);

  // project energy and water content
  double dt_next = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  double dt_prev = S_->get_time(tag_current_) - time_prev2_;

  // -- get wc and energy data
  const Epetra_MultiVector& wc0 = *wc_prev2_;
  const Epetra_MultiVector& wc1 =
    *S_->GetPtr<CompositeVector>(wc_key_, tag_current_)->ViewComponent("cell", false);
  Epetra_MultiVector& wc2 =
    *S_->GetPtrW<CompositeVector>(wc_key_, tag_next_, wc_key_)->ViewComponent("cell", false);

  const Epetra_MultiVector& e0 = *e_prev2_;
  const Epetra_MultiVector& e1 =
    *S_->GetPtr<CompositeVector>(e_key_, tag_current_)->ViewComponent("cell", false);
  Epetra_MultiVector& e2 =
    *S_->GetPtrW<CompositeVector>(e_key_, tag_next_, e_key_)->ViewComponent("cell", false);

  // -- project
  wc2 = wc0;
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  wc2.Update(dt_ratio, wc1, 1. - dt_ratio);
  e2.Update(dt_ratio, e1, 1. - dt_ratio);

  // -- extra data
  const Epetra_MultiVector& cv =
    *S_->GetPtr<CompositeVector>(cv_key_, tag_next_)->ViewComponent("cell", false);

  int rank = mesh_->getComm()->MyPID();
  int ncells = wc0.MyLength();
  for (int c = 0; c != ncells; ++c) {
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    AmanziGeometry::Point result(2);
    int ierr = 0;

    double T_guess = temp_guess_c[0][c];
    double T_prev = T1[0][c];
    double T_prev2 = (T_guess - dt_ratio * T_prev) / (1. - dt_ratio);

    double p_guess = pres_guess_c[0][c];
    double p_prev = p1[0][c];
    double p_prev2 = (p_guess - dt_ratio * p_prev) / (1. - dt_ratio);

    double p = p1[0][c];
    double T = T1[0][c];

    double wc_tmp(0.), e_tmp(0.);
    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << std::setprecision(14) << "Predicting: c = " << c << std::endl
                  << "   based upon h_old = " << dt_prev << ", h_next = " << dt_next << std::endl
                  << "   -------------" << std::endl
                  << "   Prev wc,e: " << wc1[0][c] << ", " << e1[0][c] << std::endl
                  << "   Prev p,T: " << p << ", " << T << std::endl
                  << "   -------------" << std::endl
                  << "   Extrap wc,e: " << wc2[0][c] << ", " << e2[0][c] << std::endl
                  << "   Extrap p,T: " << pres_guess_c[0][c] << ", " << T_guess << std::endl;

    model_->UpdateModel(S_.ptr(), c);
    ierr = model_->Evaluate(T_guess, pres_guess_c[0][c], e_tmp, wc_tmp);
    AMANZI_ASSERT(!ierr);
    bool ewc_completed = false;

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "   Calc wc,e of extrap: " << wc_tmp * cv[0][c] << ", " << e_tmp * cv[0][c]
                  << std::endl
                  << "   -------------" << std::endl;

    // FREEZE-THAW transition
    if (T_guess - T < 0.) { // decreasing, freezing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   decreasing temps..." << std::endl;

      if (!model_->Freezing(T_guess + cusp_size_T_freezing_, p_guess)) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
          *dcvo->os() << "   above freezing, keep T,p projections" << std::endl;
          double sl, si, sg;
          model_->EvaluateSaturations(T_guess + cusp_size_T_freezing_, p_guess, sg, sl, si);
          *dcvo->os() << "   si,sl,sg = " << si << "," << sl << "," << sg << std::endl;
        }

        // pass, guesses are good

      } else if (model_->Freezing(T_prev + cusp_size_T_freezing_, p_prev)) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   second point past the freezing point, keep T,p projections"
                      << std::endl;
        // pass, guesses are good

      } else {
        // -- invert for T,p at the projected ewc
        ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
        ewc_completed = true;
        if (ierr) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;
          // pass, keep the T,p projections
        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "     kept within the transition zone." << std::endl
                        << "   p,T = " << p << ", " << T << std::endl;

          // in the transition zone of latent heat exchange
          if (T > 200.) {
            temp_guess_c[0][c] = T;
            pres_guess_c[0][c] = p;
          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "       not admissible!" << std::endl;
          }
        }
      }
#if EWC_THAWING
    } else { // increasing, thawing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   increasing temps..." << std::endl;

      if (!model_->Freezing(T_prev + cusp_size_T_thawing_, p_prev)) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   above freezing, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else if (model_->Freezing(T_guess, p_guess)) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   projection below freezing, keep T,p projections" << std::endl;
        // pass, guesses are good

      } else {
        // in the transition zone of latent heat exchange
        ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
        ewc_completed = true;
        if (ierr) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;
          // pass, keep the T,p projections
        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "     kept within the transition zone." << std::endl
                        << "   p,T = " << p << ", " << T << std::endl;

          // two ways to get a projected T past freezing point:
          //  -- be on the lower branch and overshoot (ewc results in smaller dT)
          //  -- be on the middle branch and get over the hump (ewc results in much larger dT)
          if (T - T_prev < temp_guess_c[0][c] - T_prev) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     dT_ewc < dT_std, on the lower branch, using EWC" << std::endl;
            temp_guess_c[0][c] = T;
            pres_guess_c[0][c] = p;
          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     dT_ewc > dT_std, on the middle branch, use std prediction"
                          << std::endl;
          }
        }
      }
#endif
    }

#if EWC_SATURATION
    // SATURATED-UNSATURATED TRANSITION
    if (!ewc_completed) {     // do not do this if we already are doing ewc for temperature reasons
      if (p_guess - p < 0.) { // decreasing, becoming unsaturated
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   decreasing pressures..." << std::endl;

        if (p_guess + 100. > p_atm) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   still saturated, keep T,p projections" << std::endl;
          // pass, guesses are good

        } else if (p_prev + 100. < p_atm) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   second point past the saturation point, keep T,p projections"
                        << std::endl;
          // pass, guesses are good

        } else {
          // -- invert for T,p at the projected ewc
          ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
          if (ierr) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;
            // pass, keep the T,p projections
          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     kept within the transition zone." << std::endl
                          << "   p,T = " << p << ", " << T << std::endl;

            // in the transition zone of latent heat exchange
            if (T > 200.) {
              temp_guess_c[0][c] = T;
              pres_guess_c[0][c] = p;
            } else {
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "       not admissible!" << std::endl;
            }
          }
        }

#  if EWC_INCREASING_PRESSURE
      } else { // increasing, becoming saturated
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   increasing pressures..." << std::endl;

        if (p_prev + 10. > p_atm) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   saturated, keep T,p projections" << std::endl;
          // pass, guesses are good
        } else if (p_guess < p_atm) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   projection unsaturated, keep T,p projections" << std::endl;
          // pass, guesses are good

        } else {
          // in the transition zone of latent heat exchange
          ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
          if (ierr) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;
            // pass, keep the T,p projections
          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     kept within the transition zone." << std::endl
                          << "   p,T = " << p << ", " << T << std::endl;

            // two ways to get a projected p to saturated:
            //  -- be on the lower branch and overshoot (ewc results in smaller dp)
            //  -- be on the middle branch and get over the hump (ewc results in much larger dp)
            if (p - p_prev < pres_guess_c[0][c] - p_prev) {
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     dp_ewc < dp_std, on the lower branch, using EWC" << std::endl;
              temp_guess_c[0][c] = T;
              pres_guess_c[0][c] = p;
            } else {
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     dp_ewc > dp_std, on the middle branch, use std prediction"
                            << std::endl;
            }
          }
        }
#  endif
      }
    }
#endif
  }
  return true;
}


void
MPCDelegateEWCSubsurface::precon_ewc_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  // Pu currently stores (dp_std,dT_std).  The approach here is:
  // -- calculate (wc_ewc,e_ewc), a wc and energy-basis correction:
  //        wc_ewc = wc_0 - [ dwc/dp  dwc/dT ]( dp_std )
  //         e_ewc =  e_0 - [  de/dp   de/dT ]( dT_std )
  // -- invert the relation, calculating (p_ewc, T_ewc) = ewc^-1(wc_ewc, e_ewc)
  // -- use the smaller update of ||(p_ewc,T_ewc) - (p_0,T_0)|| and ||(dp_std, dT_std)||.

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Preconditioning using SmartEWC algorithm" << std::endl;

  // projected guesses for T and p
  Epetra_MultiVector& dp_std = *Pu->SubVector(0)->Data()->ViewComponent("cell", false);
  Epetra_MultiVector& dT_std = *Pu->SubVector(1)->Data()->ViewComponent("cell", false);

  // additional data required
  const Epetra_MultiVector& cv =
    *S_->GetPtr<CompositeVector>(cv_key_, tag_next_)->ViewComponent("cell", false);
  const double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  // old values
  const Epetra_MultiVector& p_old =
    *S_->Get<CompositeVector>(pres_key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& T_old =
    *S_->Get<CompositeVector>(temp_key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& wc_old =
    *S_->Get<CompositeVector>(wc_key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& e_old =
    *S_->Get<CompositeVector>(e_key_, tag_next_).ViewComponent("cell", false);

  // min change values... ewc is not useful near convergence
  double dT_min = 0.01;
  double dp_min = 100.;

  int rank = mesh_->getComm()->MyPID();
  int ncells = cv.MyLength();
  for (int c = 0; c != ncells; ++c) {
    // debugger
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    double T_prev = T_old[0][c];
    double T_std = T_prev - dT_std[0][c];
    double p_prev = p_old[0][c];
    double p_std = p_prev - dp_std[0][c];
    bool precon_ewc = false;

    model_->UpdateModel(S_.ptr(), c);
    bool ewc_completed = false;

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Precon: c = " << c << std::endl;

    // only do EWC if corrections are large, ie not clearly converging
    if (std::abs(dT_std[0][c]) > dT_min || std::abs(dp_std[0][c]) > dp_min) {
      // FREEZE-THAW TRANSITION
      if (-dT_std[0][c] < 0.) { // decreasing, freezing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   decreasing temps..." << std::endl;

        if (!model_->Freezing(T_std, p_std)) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   above freezing, keep std correction" << std::endl;
          // pass, guesses are good

        } else if (model_->Freezing(T_prev, p_prev)) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   linearization point past the freezing point, keep std correction"
                        << std::endl;
          // pass, guesses are good

        } else {
          // calculate the correction in ewc
          double wc_ewc =
            wc_old[0][c] - (jac_[c](0, 0) * dp_std[0][c] + jac_[c](0, 1) * dT_std[0][c]);
          double e_ewc =
            e_old[0][c] - (jac_[c](1, 0) * dp_std[0][c] + jac_[c](1, 1) * dT_std[0][c]);
          bool verbose = false;
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
            verbose = true;
            *dcvo->os() << std::setprecision(14) << "   Prev p,T: " << p_old[0][c] << ", "
                        << T_old[0][c] << std::endl
                        << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                        << "   -------------" << std::endl
                        << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c]
                        << std::endl
                        << "   applying the EWC precon:" << std::endl
                        << "     Jac = [" << jac_[c](0, 0) << ", " << jac_[c](0, 1)
                        << "] = " << wc_old[0][c] - wc_ewc << std::endl
                        << "           [" << jac_[c](1, 0) << ", " << jac_[c](1, 1)
                        << "] = " << e_old[0][c] - e_ewc << std::endl
                        << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;
          }

          // -- invert for T,p at the projected ewc
          double T(T_prev), p(p_old[0][c]);
          int ierr = model_->InverseEvaluate(e_ewc / cv[0][c], wc_ewc / cv[0][c], T, p, verbose);
          ewc_completed = true;
          if (ierr) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PRECON" << std::endl;
            // pass, keep the T,p projections

          } else {
            double dT_ewc = T_old[0][c] - T;
            double dp_ewc = p_old[0][c] - p;

            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
              *dcvo->os() << "     kept within the transition zone." << std::endl
                          << "   p,T_ewc = " << p << ", " << T << std::endl
                          << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                *dcvo->os() << "  sufficient change" << std::endl;
              } else {
                *dcvo->os() << "  insufficient change, trying anyway" << std::endl;
              }
            }

            //            if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
            dT_std[0][c] = dT_ewc;
            dp_std[0][c] = dp_ewc;
            //            }
          }
        }

#if EWC_PC_THAWING
      } else { // increasing, thawing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   increasing temps..." << std::endl;

        if (!model_->Freezing(T_prev, p_prev)) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   above freezing, keep std correction" << std::endl;
          // pass, update is are good

        } else if (model_->Freezing(T_std, p_std)) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   still frozen, keep std correction" << std::endl;
          // pass, update is are good

        } else {
          // calculate the correction in ewc
          double wc_ewc =
            wc_old[0][c] - (jac_[c](0, 0) * dp_std[0][c] + jac_[c](0, 1) * dT_std[0][c]);
          double e_ewc =
            e_old[0][c] - (jac_[c](1, 0) * dp_std[0][c] + jac_[c](1, 1) * dT_std[0][c]);
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   Prev p,T: " << p_old[0][c] << ", " << T_old[0][c] << std::endl
                        << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                        << "   -------------" << std::endl
                        << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c]
                        << std::endl
                        << "   applying the EWC precon:" << std::endl
                        << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;

          // -- invert for T,p at the projected ewc
          double T(T_prev), p(p_old[0][c]);
          int ierr = model_->InverseEvaluate(e_ewc / cv[0][c], wc_ewc / cv[0][c], T, p);
          ewc_completed = true;

          if (ierr) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PRECON" << std::endl;
            // pass, keep the T,p projections

          } else {
            double dT_ewc = T_old[0][c] - T;
            double dp_ewc = p_old[0][c] - p;

            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
              *dcvo->os() << "     within the transition zone." << std::endl
                          << "   p,T_ewc = " << p << ", " << T << std::endl
                          << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
              if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) &&
                  (std::abs(dT_ewc) < std::abs(dT_std[0][c]))) {
                *dcvo->os()
                  << "  sufficient change, and decreased dT (and so on the lower branch), using EWC"
                  << std::endl;
              } else if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min)) {
                *dcvo->os() << "  increased dT (and so on the middle branch), using std"
                            << std::endl;
              } else {
                *dcvo->os() << "  insufficient change, trying anyway" << std::endl;
              }
            }

            //            if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) &&
            if (std::abs(dT_ewc) < std::abs(dT_std[0][c])) {
              dT_std[0][c] = dT_ewc;
              dp_std[0][c] = dp_ewc;
            }
          }
        }
#endif
      }

#if EWC_PC_SATURATION
      if (!ewc_completed) {
        // SATURATED-UNSATURATED TRANSITION
        if (-dp_std[0][c] < 0.) { // decreasing, going unsaturated
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   decreasing pressures..." << std::endl;

          if (p_std + 100. > p_atm) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   guess saturated, keep std correction" << std::endl;
            // pass, guesses are good

          } else if (p_prev + 100. < p_atm) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   linearization point past p_atm, keep std correction" << std::endl;
            // pass, guesses are good

          } else {
            // calculate the correction in ewc
            double wc_ewc =
              wc_old[0][c] - (jac_[c](0, 0) * dp_std[0][c] + jac_[c](0, 1) * dT_std[0][c]);
            double e_ewc =
              e_old[0][c] - (jac_[c](1, 0) * dp_std[0][c] + jac_[c](1, 1) * dT_std[0][c]);
            bool verbose = false;
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
              verbose = true;
              *dcvo->os() << std::setprecision(14) << "   Prev p,T: " << p_old[0][c] << ", "
                          << T_old[0][c] << std::endl
                          << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                          << "   -------------" << std::endl
                          << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c]
                          << std::endl
                          << "   applying the EWC precon:" << std::endl
                          << "     Jac = [" << jac_[c](0, 0) << ", " << jac_[c](0, 1)
                          << "] = " << wc_old[0][c] - wc_ewc << std::endl
                          << "           [" << jac_[c](1, 0) << ", " << jac_[c](1, 1)
                          << "] = " << e_old[0][c] - e_ewc << std::endl
                          << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;
            }

            // -- invert for T,p at the projected ewc
            double T(T_prev), p(p_old[0][c]);
            int ierr = model_->InverseEvaluate(e_ewc / cv[0][c], wc_ewc / cv[0][c], T, p, verbose);
            ewc_completed = true;
            if (ierr) {
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "FAILED EWC PRECON" << std::endl;
              // pass, keep the T,p projections

            } else {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     kept within the transition zone." << std::endl
                            << "   p,T_ewc = " << p << ", " << T << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                  *dcvo->os() << "  sufficient change" << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change, taking anyway" << std::endl;
                }
              }

              //              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
              dT_std[0][c] = dT_ewc;
              dp_std[0][c] = dp_ewc;
              //              }
            }
          }

#  if EWC_PC_INCREASING_PRESSURE
        } else { // increasing, thawing
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   increasing pressures..." << std::endl;

          if (p_prev > p_atm) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   saturated, keep std correction" << std::endl;
            // pass, update is are good

          } else if (p_std < p_atm) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   still unsaturated, keep std correction" << std::endl;
            // pass, update is are good

          } else {
            // calculate the correction in ewc
            double wc_ewc =
              wc_old[0][c] - (jac_[c](0, 0) * dp_std[0][c] + jac_[c](0, 1) * dT_std[0][c]);
            double e_ewc =
              e_old[0][c] - (jac_[c](1, 0) * dp_std[0][c] + jac_[c](1, 1) * dT_std[0][c]);
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   Prev p,T: " << p_old[0][c] << ", " << T_old[0][c] << std::endl
                          << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                          << "   -------------" << std::endl
                          << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c]
                          << std::endl
                          << "   applying the EWC precon:" << std::endl
                          << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;

            // -- invert for T,p at the projected ewc
            double T(T_prev), p(p_old[0][c]);
            int ierr = model_->InverseEvaluate(e_ewc / cv[0][c], wc_ewc / cv[0][c], T, p);

            if (ierr) {
              ewc_completed = true;
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "FAILED EWC PRECON" << std::endl;
              // pass, keep the T,p projections

            } else {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     within the transition zone." << std::endl
                            << "   p,T_ewc = " << p << ", " << T << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) &&
                    (std::abs(dp_ewc) < std::abs(dp_std[0][c]))) {
                  *dcvo->os() << "  sufficient change, and decreased dp (and so on the lower "
                                 "branch), using EWC"
                              << std::endl;
                } else if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min)) {
                  *dcvo->os() << "  increased dp (and so on the middle branch), using std"
                              << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change, taking anyway" << std::endl;
                }
              }

              //              if ((std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) &&
              if (std::abs(dp_ewc) < std::abs(dp_std[0][c])) {
                ewc_completed = true;
                dT_std[0][c] = dT_ewc;
                dp_std[0][c] = dp_ewc;
              }
            }
          }
#  endif
        }
      }
#endif
    }
  }
}

} // namespace Amanzi
