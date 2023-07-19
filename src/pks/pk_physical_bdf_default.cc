/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Standard base for most PKs, this combines both domains/meshes of
PKPhysicalBase and BDF methods of PK_BDF_Default.
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"
#include "pk_helpers.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void
PK_PhysicalBDF_Default::Setup()
{
  // call the meat of the base constructurs via Setup methods
  PK_Physical_Default::Setup();
  PK_BDF_Default::Setup();

  // boundary conditions
  bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));

  // convergence criteria is based on a conserved quantity
  if (conserved_key_.empty()) {
    conserved_key_ = Keys::readKey(*plist_, domain_, "conserved quantity");
  }
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, true);
  // we also use a copy of the conserved quantity, as this is a better choice in the error norm
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_, true);

  // cell volume used throughout
  if (cell_vol_key_.empty()) {
    cell_vol_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  }
  requireAtNext(cell_vol_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, true);

  atol_ = plist_->get<double>("absolute error tolerance", 1.0);
  rtol_ = plist_->get<double>("relative error tolerance", 1.0);
  fluxtol_ = plist_->get<double>("flux error tolerance", 1.0);
};


// -----------------------------------------------------------------------------
// initialize.  Note both BDFBase and PhysicalBase have initialize()
// methods, so we need a unique overrider.
// -----------------------------------------------------------------------------
void
PK_PhysicalBDF_Default::Initialize()
{
  // Just calls both subclass's initialize.  NOTE - order is important here --
  // PhysicalBase grabs the primary variable and stuffs it into the solution,
  // which must be done prior to BDFBase initializing the timestepper.
  PK_Physical_Default::Initialize();
  PK_BDF_Default::Initialize();
}


int
PK_PhysicalBDF_Default::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                            Teuchos::RCP<TreeVector> Pu)
{
  *Pu = *u;
  return 0;
}


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double
PK_PhysicalBDF_Default::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<const TreeVector> res)
{
  // Abs tol based on old conserved quantity -- we know these have been vetted
  // at some level whereas the new quantity is some iterate, and may be
  // anything from negative to overflow.
  //  S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name()); // for the future...
  const Epetra_MultiVector& conserved =
    *S_->Get<CompositeVector>(conserved_key_, tag_current_).ViewComponent("cell", true);
  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = res->Data();
  double h = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp = dvec->begin(); comp != dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == "cell") {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp, false);
      for (unsigned int c = 0; c != ncells; ++c) {
        double enorm_c =
          std::abs(h * dvec_v[0][c]) / (atol_ * cv[0][c] + rtol_ * std::abs(conserved[0][c]));
        AMANZI_ASSERT((atol_ * cv[0][c] + rtol_ * std::abs(conserved[0][c])) > 0.);

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }

    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f = 0; f != nfaces; ++f) {
        auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::OWNED);
        double cv_min =
          cells.size() == 1 ? cv[0][cells[0]] : std::min(cv[0][cells[0]], cv[0][cells[1]]);
        double conserved_min = cells.size() == 1 ?
                                 conserved[0][cells[0]] :
                                 std::min(conserved[0][cells[0]], conserved[0][cells[1]]);

        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f]) /
                         (atol_ * cv_min + rtol_ * std::abs(conserved_min));
        AMANZI_ASSERT((atol_ * cv_min + rtol_ * std::abs(conserved_min)) > 0.);
        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      // double norm;
      // dvec_v.Norm2(&norm);
      //      AMANZI_ASSERT(norm < 1.e-15);
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr;

      ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      AMANZI_ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] ("
                 << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;

  int ierr;
  ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, comm);
  AMANZI_ASSERT(!ierr);
  return enorm_val;
};


void
PK_PhysicalBDF_Default::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  PK_BDF_Default::CommitStep(t_old, t_new, tag_next);
  PK_Physical_Default::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  // copy over conserved quantity
  assign(conserved_key_, tag_current, tag_next, *S_);
}


void
PK_PhysicalBDF_Default::FailStep(double t_old, double t_new, const Tag& tag)
{
  PK_Physical_Default::FailStep(t_old, t_new, tag);
}


// -----------------------------------------------------------------------------
// Calling this indicates that the time integration scheme is changing the
// value of the solution in state.
// -----------------------------------------------------------------------------
void
PK_PhysicalBDF_Default::ChangedSolution(const Tag& tag)
{
  Teuchos::RCP<Evaluator> fm = S_->GetEvaluatorPtr(key_, tag);
  Teuchos::RCP<EvaluatorPrimaryCV> solution_evaluator =
    Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(fm);
  AMANZI_ASSERT(solution_evaluator != Teuchos::null);
  solution_evaluator->SetChanged();
};

} // namespace Amanzi
