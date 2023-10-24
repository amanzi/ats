/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Implementation for the NullEnergy PK.  This PK simply provides a constant
   temperature, and is provided for testing with other PKs that depend upon an
   energy equation.  This could easily be provided by the state as an independent
   variable, but this is nice for testing the full hierarchy with a simple PK.

   Example usage:

   <ParameterList name="energy">
   <Parameter name="PK model" type="string" value="Constant Temperature"/>
   <Parameter name="Constant Temperature" type="double" value="290.0"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "test_snow_dist.hh"

namespace Amanzi {

void
TestSnowDist::setup(const Teuchos::Ptr<State>& S)
{
  PKPhysicalBase::setup(S);

  S->Require<CompositeVector, CompositeVectorSpace>(key_, Tags::NEXT, name_)
    .SetMesh(S->GetMesh("surface"))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S->RequireEvaluator("precipitation_snow");
  S->Require<CompositeVector, CompositeVectorSpace>("precipitation_snow", Tags::NEXT)
    .SetMesh(S->GetMesh("surface"))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
};

// -- call your favorite
bool
TestSnowDist::advance(double dt)
{
  if (sink_type_ == "factor") {
    S_next_->GetW<CompositeVector>("snow_depth", name_).Scale(sink_value_);
  } else if (sink_type_ == "constant") {
    Epetra_MultiVector& sd =
      *S_next_->GetW<CompositeVector>("snow_depth", name_).ViewComponent("cell", false);
    for (int c = 0; c != sd.MyLength(); ++c) {
      sd[0][c] -= 10 * dt * sink_value_;
      sd[0][c] = std::max(0., sd[0][c]);
    }
  }

  S_next_->GetEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetPtrW<CompositeVector>("snow_depth", name_)
    ->Update(10. * dt,
             *S_next_->GetPtr<CompositeVector>("precipitation_snow"),
             1.); // factor of 10 for SWE-to-snow ht conversion

  solution_evaluator_->SetChanged(S_next_.ptr());
  return false;
};


} // namespace Amanzi
