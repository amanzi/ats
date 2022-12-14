/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

{
  if_elseif
}
for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
  {
    { keyEpetraVectorList } Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

    int ncomp = result->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      {
        result_v[0][i] = model_->D{ myKeyMethod } D{ wrtMethod }({ myMethodArgs });
      }
    }
  }
}
