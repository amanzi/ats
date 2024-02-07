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

#include "ewc_model.hh"
#include "mpc_delegate_ewc_surface.hh"

namespace Amanzi {

MPCDelegateEWCSurface::MPCDelegateEWCSurface(Teuchos::ParameterList& plist,
                                             const Teuchos::RCP<State>& S)
  : MPCDelegateEWC(plist, S)
{
  T_cutoff_ = plist_->get<double>("fully frozen temperature", 272.15);
}

bool
MPCDelegateEWCSurface::modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell", false);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "  Modifying surface predictor using SmartEWC algorithm" << std::endl;
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

  // Ensure the necessity of doing this... if max(pres_guess_c) < p_atm then there is no water anywhere.
  double p_max;
  p1.MaxValue(&p_max);
  if (p_max < 101325.) { return false; }

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

    double p = p1[0][c];
    double T = T1[0][c];

    double wc_tmp(0.), e_tmp(0.);
    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Predicting: sc = " << c << std::endl
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

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "   Calc wc,e of extrap: " << wc_tmp * cv[0][c] << ", " << e_tmp * cv[0][c]
                  << std::endl
                  << "   -------------" << std::endl;

    if (T_guess - T < 0.) { // decreasing, freezing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   decreasing temps..." << std::endl;

      if (T_guess >= 273.149) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   above freezing range, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else if (T_prev2 < T_cutoff_) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   below freezing range, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else {
        // -- invert for T,p at the projected ewc
        ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
        if (ierr) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;

        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   new p,T projection: " << p << ", " << T << std::endl;

          // in the transition zone of latent heat exchange
          temp_guess_c[0][c] = T;
          pres_guess_c[0][c] = p;
        }
      }

    } else { // increasing, thawing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   increasing temps..." << std::endl;

      if (T_prev2 >= 273.149) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   above freezing range, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else if (T_guess < T_cutoff_) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   below freezing range, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else {
        // in the transition zone of latent heat exchange
        ierr = model_->InverseEvaluate(e2[0][c] / cv[0][c], wc2[0][c] / cv[0][c], T, p);
        if (ierr) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;
        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   new p,T projection: " << p << ", " << T << std::endl;

          temp_guess_c[0][c] = T;
          pres_guess_c[0][c] = p;
        }
      }
    }
  }
  return true;
}


void
MPCDelegateEWCSurface::precon_ewc_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
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
    *S_->Get<CompositeVector>("cell_volume", tag_next_).ViewComponent("cell", false);

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
    bool precon_ewc = false;

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Precon: sc = " << c << std::endl;

    // only do EWC if corrections are large, ie not clearly converging
    model_->UpdateModel(S_.ptr(), c);
    if (std::abs(dT_std[0][c]) > dT_min || std::abs(dp_std[0][c]) > dp_min) {
      if (-dT_std[0][c] < 0.) { // decreasing, freezing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   decreasing temps..." << std::endl;

        if (T_std >= 273.15) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   above freezing, keep std correction" << std::endl;
          // pass, guesses are good

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
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PRECON" << std::endl;

          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   new p,T projection: " << p << ", " << T << std::endl;

            // Check if our energy projection is before the 2nd inflection point
            if (T > T_cutoff_) {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     kept within the transition zone." << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                  *dcvo->os() << "  sufficient change" << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change" << std::endl;
                }
              }

              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                dT_std[0][c] = dT_ewc;
                dp_std[0][c] = dp_ewc;
              }
            } else {
              // pass, past the transition and we are just chilling ice
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     outside the transition zone." << std::endl;
            }
          }
        }
      } else { // increasing, thawing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   increasing temps..." << std::endl;

        if (T_std >= 273.15) {
          // pass, update is are good
        } else {
          if (T_std < T_cutoff_) {
            // pass, warming ice but still outside of the transition region
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   above freezing, keep T,p projections" << std::endl;
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

            if (!ierr && T < 273.15) {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     kept within the transition zone." << std::endl
                            << "   p,T_ewc = " << p << ", " << T << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                  *dcvo->os() << "  sufficient change" << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change" << std::endl;
                }
              }

              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                dT_std[0][c] = dT_ewc;
                dp_std[0][c] = dp_ewc;
              }

            } else {
              // pass, past the transition zone and onto the upper branch
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     outside the transition zone." << std::endl;
            }
          }
        }
      }
    }
  }
}

} // namespace Amanzi
