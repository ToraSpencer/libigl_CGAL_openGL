#ifndef IGL_DOUBLEAREA_H
#define IGL_DOUBLEAREA_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // DOUBLEAREA computes twice the area for each input triangle[quad]
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for vers (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for tris (e.g. derived from
  //     MatrixXi)
  //   DeriveddblA  derived type of eigen matrix for dblA (e.g. derived from
  //     MatrixXd)
  // Inputs:
  //   vers  #vers by dim list of mesh vertex positions
  //   tris  #tris by simplex_size list of mesh faces (must be triangles or quads)
  // Outputs:
  //   dblA  #tris list of triangle[quad] double areas (SIGNED only for 2D input)
  //
  // Known bug: For dim==3 complexity is O(#vers + #tris)!! Not just O(#tris). This is a big deal
  // if you have 1million unreferenced vertices and 1 face
  template <typename DerivedV, typename DerivedF, typename DeriveddblA>
  IGL_INLINE void doublearea(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DeriveddblA> & dblA);
  // Stream of triangles, computes signed area...
  template <
    typename DerivedA,
    typename DerivedB,
    typename DerivedC,
    typename DerivedD>
  IGL_INLINE void doublearea(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B,
    const Eigen::MatrixBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedD> & D);
  // Single triangle in 2D!
  //
  // This should handle streams of corners not just single corners
  template <
    typename DerivedA,
    typename DerivedB,
    typename DerivedC>
  IGL_INLINE typename DerivedA::Scalar doublearea_single(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B,
    const Eigen::MatrixBase<DerivedC> & C);
  // Same as above but use instrinsic edge lengths rather than (vers,tris) mesh. This
  //
  // Inputs:
  //   l  #tris by dim list of edge lengths using
  //     for triangles, columns correspond to edges 23,31,12
  //   nan_replacement  what value should be used for triangles whose given
  //     edge lengths do not obey the triangle inequality. These may be very
  //     wrong (e.g., [100 1 1]) or may be nearly degenerate triangles whose
  //     floating point side length computation leads to breach of the triangle
  //     inequality. One may wish to set this parameter to 0 if side lengths l
  //     are _known_ to come from a valid embedding (e.g., some mesh (vers,tris)). In
  //     that case, the only circumstance the triangle inequality is broken is
  //     when the triangle is nearly degenerate and floating point error
  //     dominates: hence replacing with zero is reasonable.
  // Outputs:
  //   dblA  #tris list of triangle double areas
  template <typename Derivedl, typename DeriveddblA>
  IGL_INLINE void doublearea(
    const Eigen::MatrixBase<Derivedl> & l,
    const typename Derivedl::Scalar nan_replacement,
    Eigen::PlainObjectBase<DeriveddblA> & dblA);
  // default behavior is to assert on NaNs and leave them in place
  template <typename Derivedl, typename DeriveddblA>
  IGL_INLINE void doublearea(
    const Eigen::MatrixBase<Derivedl> & l,
    Eigen::PlainObjectBase<DeriveddblA> & dblA);
  // DOUBLEAREA_QUAD computes twice the area for each input quadrilateral
  //
  // Inputs:
  //   vers  #vers by dim list of mesh vertex positions
  //   tris  #tris by simplex_size list of mesh faces (must be quadrilaterals)
  // Outputs:
  //   dblA  #tris list of quadrilateral double areas
  //
  template <typename DerivedV, typename DerivedF, typename DeriveddblA>
  IGL_INLINE void doublearea_quad(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  Eigen::PlainObjectBase<DeriveddblA> & dblA);

}

#ifndef IGL_STATIC_LIBRARY
#  include "doublearea.cpp"
#endif

#endif
