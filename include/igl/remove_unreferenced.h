#ifndef IGL_REMOVE_UNREFERENCED_H
#define IGL_REMOVE_UNREFERENCED_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // Remove unreferenced vertices from V, updating F accordingly
  //
  // Input:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by ss list of simplices (Values of -1 are quitely skipped)
  // Outputs:
  //   NV  #NV by dim list of mesh vertex positions
  //   NF  #NF by ss list of simplices
  //   I   #V by 1 list of indices such that: NF = IM(F) and NT = IM(T)
  //      and V(find(IM<=size(NV,1)),:) = NV
  //   J  #NV by 1 list, such that NV = V(J,:)
  //   
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI>
  IGL_INLINE void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedNV> &NV,
    Eigen::PlainObjectBase<DerivedNF> &NF,
    Eigen::PlainObjectBase<DerivedI> &I);
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI,
    typename DerivedJ>
  IGL_INLINE void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedNV> &NV,
    Eigen::PlainObjectBase<DerivedNF> &NF,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);
  // Inputs:
  //   n  number of vertices (possibly greater than F.maxCoeff()+1)
  //   F  #F by ss list of simplices
  // Outputs:
  //   IM  #V by 1 list of indices such that: NF = IM(F) and NT = IM(T)
  //      and V(find(IM<=size(NV,1)),:) = NV
  //   J  #RV by 1 list, such that RV = V(J,:)
  //   
  template <
    typename DerivedF,
    typename DerivedI,
    typename DerivedJ>
  IGL_INLINE void remove_unreferenced(
    const size_t n,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);

}

#ifndef IGL_STATIC_LIBRARY
#  include "remove_unreferenced.cpp"
#endif

#endif
