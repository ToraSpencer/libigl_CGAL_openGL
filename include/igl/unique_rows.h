#ifndef IGL_UNIQUE_ROWS_H
#define IGL_UNIQUE_ROWS_H
#include "igl_inline.h"

#include <vector>
#include <Eigen/Core>

namespace igl
{
  // Act like matlab's [C,IA,IC] = unique(X,'rows')
  /*
       Templates:
         DerivedA derived scalar type, e.g. MatrixXi or MatrixXd
         DerivedIA derived integer type, e.g. MatrixXi
         DerivedIC derived integer type, e.g. MatrixXi
       Inputs:
         A  m by n matrix whose entries are to unique'd according to rows
       Outputs:
         C  #C vector of unique rows in A
         IA  #C index vector so that C = A(IA,:);
         IC  #A index vector so that A = C(IC,:);
   */
  template <typename DerivedA, typename DerivedC, typename DerivedIA, typename DerivedIC>
  IGL_INLINE void unique_rows(
    const Eigen::DenseBase<DerivedA>& A,
    Eigen::PlainObjectBase<DerivedC>& C,
    Eigen::PlainObjectBase<DerivedIA>& IA,
    Eigen::PlainObjectBase<DerivedIC>& IC);

}

#ifndef IGL_STATIC_LIBRARY
#  include "unique_rows.cpp"
#endif

#endif
