#ifndef IGL_CONNECT_BOUNDARY_TO_INFINITY_H
#define IGL_CONNECT_BOUNDARY_TO_INFINITY_H
#include "igl_inline.h"
#include <Eigen/Core>


namespace igl
{
  // Connect all boundary edges to a fictitious point at infinity.
  //
  // Inputs:
  //   tris  #tris by 3 list of face indices into some vers
  // Outputs:
  //   FO  #tris+#O by 3 list of face indices into [vers;inf inf inf], original tris are
  //     guaranteed to come first. If (vers,tris) was a manifold mesh, now it is
  //     closed with a possibly non-manifold vertex at infinity (but it will be edge-manifold).
  template <typename DerivedF, typename DerivedFO>
  IGL_INLINE void connect_boundary_to_infinity(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DerivedFO> & FO);


  // Inputs:
  //   inf_index  index of point at infinity (usually vers.rows() or tris.maxCoeff())
  template <typename DerivedF, typename DerivedFO>
  IGL_INLINE void connect_boundary_to_infinity(
    const Eigen::MatrixBase<DerivedF> & tris,
    const typename DerivedF::Scalar inf_index,
    Eigen::PlainObjectBase<DerivedFO> & FO);
  // Inputs:
  //   vers  #vers by 3 list of vertex positions
  //   tris  #tris by 3 list of face indices into some vers
  // Outputs:
  //   VO  #vers+1 by 3 list of vertex positions, original vers are guaranteed to
  //     come first. Last point is inf, inf, inf
  //   FO  #tris+#O by 3 list of face indices into VO
   

  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedVO, 
    typename DerivedFO>
  IGL_INLINE void connect_boundary_to_infinity(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DerivedVO> & VO,
    Eigen::PlainObjectBase<DerivedFO> & FO);
}

#ifndef IGL_STATIC_LIBRARY
#  include "connect_boundary_to_infinity.cpp"
#endif
#endif
