#ifndef IGL_REMOVE_UNREFERENCED_H
#define IGL_REMOVE_UNREFERENCED_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // Remove unreferenced vertices from vers, updating tris accordingly
  //
  // Input:
  //   vers  #vers by dim list of mesh vertex positions
  //   tris  #tris by ss list of simplices (Values of -1 are quitely skipped)
  // Outputs:
  //   NV  #NV by dim list of mesh vertex positions
  //   NF  #NF by ss list of simplices
  //   I   #vers by 1 list of indices such that: NF = IM(tris) and NT = IM(T)
  //      and vers(find(IM<=size(NV,1)),:) = NV
  //   J  #NV by 1 list, such that NV = vers(J,:)
  //   
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI>
  IGL_INLINE void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &vers,
    const Eigen::MatrixBase<DerivedF> &tris,
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
    const Eigen::MatrixBase<DerivedV> &vers,
    const Eigen::MatrixBase<DerivedF> &tris,
    Eigen::PlainObjectBase<DerivedNV> &NV,
    Eigen::PlainObjectBase<DerivedNF> &NF,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);
  // Inputs:
  //   n  number of vertices (possibly greater than tris.maxCoeff()+1)
  //   tris  #tris by ss list of simplices
  // Outputs:
  //   IM  #vers by 1 list of indices such that: NF = IM(tris) and NT = IM(T)
  //      and vers(find(IM<=size(NV,1)),:) = NV
  //   J  #RV by 1 list, such that RV = vers(J,:)
  //   
  template <
    typename DerivedF,
    typename DerivedI,
    typename DerivedJ>
  IGL_INLINE void remove_unreferenced(
    const size_t n,
    const Eigen::MatrixBase<DerivedF> &tris,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);

}

#ifndef IGL_STATIC_LIBRARY
#  include "remove_unreferenced.cpp"
#endif

#endif
