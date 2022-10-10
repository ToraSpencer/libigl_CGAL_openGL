#ifndef IGL_BARYCENTER_H
#define IGL_BARYCENTER_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
    // 计算每个单纯形的重心：
  /*
      Computes the barycenter of every simplex
  
       Inputs:
         vers         #vers x dim matrix of vertex coordinates
         tris         #tris x simplex_size  matrix of indices of simplex corners into vers

       Output:
         baryCenters  #tris x dim matrix of 3d vertices
  */
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedBC>
  IGL_INLINE void barycenter(
      const Eigen::MatrixBase<DerivedV> & vers,
      const Eigen::MatrixBase<DerivedF> & tris,
      Eigen::PlainObjectBase<DerivedBC> & baryCenters);
}

#ifndef IGL_STATIC_LIBRARY
#  include "barycenter.cpp"
#endif

#endif
