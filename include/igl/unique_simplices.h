#ifndef IGL_UNIQUE_SIMPLICES_H
#define IGL_UNIQUE_SIMPLICES_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
    // unique_simplices()¡ª¡ª
    /*
       Find *combinatorially* unique simplices in tris.  **Order independent**
  
       Inputs:
                 tris  #tris by simplex-size list of simplices
       Outputs:
                 FF  #FF by simplex-size list of unique simplices in tris
                 IA  #FF index vector so that FF == sort(tris(IA,:),2);
                 IC  #tris index vector so that sort(tris,2) == FF(IC,:);
    */
  template <
    typename DerivedF,
    typename DerivedFF,
    typename DerivedIA,
    typename DerivedIC>
  IGL_INLINE void unique_simplices(
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedFF>& FF,
    Eigen::PlainObjectBase<DerivedIA>& IA,
    Eigen::PlainObjectBase<DerivedIC>& IC);
  template <
    typename DerivedF,
    typename DerivedFF>
  IGL_INLINE void unique_simplices(
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedFF>& FF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "unique_simplices.cpp"
#endif

#endif
