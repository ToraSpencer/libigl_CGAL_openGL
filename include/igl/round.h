#ifndef IGL_ROUND_H
#define IGL_ROUND_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Round a scalar value
  //
  // Inputs:
  //   x  number
  // Returns x rounded to integer
  template <typename DerivedX>
  DerivedX round(const DerivedX r);


  // Round a given matrix to nearest integers
  //
  // Inputs:
  //   X  m by n matrix of scalars
  // Outputs:
  //   Y  m by n matrix of rounded integers
  template < typename DerivedX, typename DerivedY>
  IGL_INLINE void round(
    const Eigen::PlainObjectBase<DerivedX>& X,
    Eigen::PlainObjectBase<DerivedY>& Y);
}

#ifndef IGL_STATIC_LIBRARY
#  include "round.cpp"
#endif

#endif
