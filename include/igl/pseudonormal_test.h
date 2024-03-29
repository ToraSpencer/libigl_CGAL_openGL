#ifndef IGL_PSEUDONORMAL_TEST_H
#define IGL_PSEUDONORMAL_TEST_H
#include "igl_inline.h"
#include <Eigen/Core>


namespace igl
{
    /*
       Given a mesh (vers,tris), a query point q, and a point on (vers,tris) c, determine
       whether q is inside (vers,tris) --> s=-1 or outside (vers,tris) s=1, based on the
       sign of the dot product between (q-c) and normal, where normal is the normal _at c_,
       carefully chosen according to [Bærentzen & Aanæs 2005]
  
       Inputs:
         vers  #vers by 3 list of vertex positions
         tris  #tris by 3 list of triangle indices
         triNorms  #tris by 3 list of triangle normals 
         verNorms  #vers by 3 list of vertex normals (ANGLE WEIGHTING)
         edgeNorms  #E by 3 list of edge normals (UNIFORM WEIGHTING)
         edgeUeInfo  #tris*3 mapping edges in tris to E
         q          Query point
         triIdx0  index into tris to face to which c belongs
         c      Point on (vers, tris) 

       Outputs:
         s  sign
         normal  normal
  */
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedFN,
    typename DerivedVN,
    typename DerivedEN,
    typename DerivedEMAP,
    typename Derivedq,
    typename Derivedc,
    typename Scalar,
    typename Derivedn>
  IGL_INLINE void pseudonormal_test(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedFN> & triNorms,
    const Eigen::MatrixBase<DerivedVN> & verNorms,
    const Eigen::MatrixBase<DerivedEN> & edgeNorms,
    const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo,
    const Eigen::MatrixBase<Derivedq> & q,
    const int triIdx0,
    Eigen::PlainObjectBase<Derivedc> & c,
    Scalar & s,
    Eigen::PlainObjectBase<Derivedn> & normal);


  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedEN,
    typename DerivedVN,
    typename Derivedq,
    typename Derivedc,
    typename Scalar,
    typename Derivedn>
  IGL_INLINE void pseudonormal_test(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & E,
    const Eigen::MatrixBase<DerivedEN> & edgeNorms,
    const Eigen::MatrixBase<DerivedVN> & verNorms,
    const Eigen::MatrixBase<Derivedq> & q,
    const int e,
    Eigen::PlainObjectBase<Derivedc> & c,
    Scalar & s,
    Eigen::PlainObjectBase<Derivedn> & normal);
}
#ifndef IGL_STATIC_LIBRARY
#  include "pseudonormal_test.cpp"
#endif
#endif
