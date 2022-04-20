#ifndef _DBC_H_
#define _DBC_H_

#include "exceptions.hh"

namespace DBC {

/* DBC_assertion
 *
 * An exception class for DBC asserion volations.
 * 
 */
class Assertion : public Exceptions::Amanzi_exception {
public:
  Assertion (const char* condition, const char* file, unsigned int line);
  const char* what () const throw ();

 public:
  const char* assertion_;
  const char* filename_;
  unsigned int line_number_;
};

void amanzi_assert(const char * cond, const char * file, unsigned int line);

} // namespace DBC


// ---------------------------------
// Macros for DBC checks in the code
// ---------------------------------

// The do wrapper prevents the if statement from grabbing subsequent
// else statements away from enclosing ifs.  The version when DBC is
// not enabled compiles away to nothing, but surpresses warning about
// unused variables in the expression a.

#ifdef ENABLE_DBC
#define AMANZI_ASSERT(bool_expression) do { if (!(bool_expression)) DBC::amanzi_assert( #bool_expression, __FILE__, __LINE__ ); } while (0)
#else
#define AMANZI_ASSERT(a) do { (void)sizeof(a); } while (0);
#endif /* ENABLE_DBC */


#endif /* _DBC_H_ */

