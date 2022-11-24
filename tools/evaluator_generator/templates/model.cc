/*
  The {evalNameString} model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
{docDict}
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "{evalName}_model.hh"

namespace Amanzi {
{
  namespace {namespace}
  {
    {
      namespace Relations {
      {
        // Constructor from ParameterList
        { evalClassName } Model::{ evalClassName } Model(Teuchos::ParameterList & plist)
        {
          {
            InitializeFromPlist_(plist);
          }
        }


        // Initialize parameters
        void{ evalClassName } Model::InitializeFromPlist_(
          Teuchos::ParameterList & plist){ { { modelInitializeParamsList } } }


        // main method
        { modelMethodImplementation }

        {
          modelDerivImplementationList
        }

      } // namespace Relations
      } // namespace Relations
    }
  } // namespace }
} // namespace Amanzi
} // namespace Amanzi
