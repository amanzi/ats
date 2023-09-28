/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  The {evalNameString} evaluator is an algebraic evaluator of a given model.
{docDict}
  Generated via evaluator_generator.
*/

#include "{evalName}_evaluator.hh"
#include "{evalName}_model.hh"

namespace Amanzi {
{
  namespace {namespace}
  {
    {
      namespace Relations {
      {
        // Constructor from ParameterList
        { evalClassName } Evaluator::{ evalClassName } Evaluator(Teuchos::ParameterList & plist)
          : EvaluatorSecondaryMonotypeCV(plist)
        {
          {
            Teuchos::ParameterList& sublist = plist_.sublist("{evalName} parameters");
            model_ = Teuchos::rcp(new { evalClassName } Model(sublist));
            InitializeFromPlist_();
          } // namespace Relations
        }   // namespace Relations


        // Copy constructor
        { evalClassName } Evaluator::{ evalClassName } Evaluator(const { evalClassName } Evaluator &
                                                                 other)
          : EvaluatorSecondaryMonotypeCV(other), { keyCopyConstructorList } model_(other.model_)
        {
          {}
        }


        // Virtual copy constructor
        Teuchos::RCP<Evaluator>{ evalClassName } Evaluator::Clone() const
        {
          {
            return Teuchos::rcp(new { evalClassName } Evaluator(*this));
          }
        }


        // Initialize by setting up dependencies
        void{ evalClassName } Evaluator::InitializeFromPlist_()
        {
          {
            // Set up my dependencies
            // - defaults to prefixed via domain
            Key domain_name = Keys::getDomain(my_key_);

            // - pull Keys from plist
            {
              keyInitializeList
            }
          }
        }


        void{ evalClassName } Evaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                                        const Teuchos::Ptr<CompositeVector>& result)
        {
          {
            {
              keyCompositeVectorList
            }

            {
              evaluateModel
            }
          }
        }


        void{ evalClassName } Evaluator::EvaluateFieldPartialDerivative_(
          const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
        {
          {
            {
              keyCompositeVectorList
            }

            {
              evaluateDerivs
            }
          }
        }
      } // namespace Relations
      } // namespace Relations
    }   // namespace Amanzi
  }     // namespace }
} // namespace Amanzi
} // namespace Amanzi
