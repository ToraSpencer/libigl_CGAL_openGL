#ifndef IGL_ORIENTED_FACETS_H
#define IGL_ORIENTED_FACETS_H
#include "igl_inline.h"
#include <Eigen/Dense>


namespace igl
{
    // 计算网格的半边（有向边）
    /*
  // ORIENTED_FACETS Determines all "directed
  // [facets](https://en.wikipedia.org/wiki/Simplex#Elements)" of a given set of
  // simplicial elements. For a manifold triangle mesh, this computes all
  // half-edges. For a manifold tetrahedral mesh, this computes all half-faces.
  //
  // Inputs:
  //   F  #F by simplex_size  list of simplices
  // Outputs:
  //   E  #E by simplex_size-1  list of facets, such that E.row(f+#F*c) is the
  //     facet opposite F(f,c)
  //
  // Note: this is not the same as igl::edges because this includes every
  // directed edge including repeats (meaning interior edges on a surface will
  // show up once for each direction and non-manifold edges may appear more than
  // once for each direction).
  //
  // Note: This replaces the deprecated `all_edges` function
    */
  template <typename DerivedF, typename DerivedE>
  IGL_INLINE void oriented_facets(
    const Eigen::MatrixBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedE> & E);

}

#ifndef IGL_STATIC_LIBRARY
#  include "oriented_facets.cpp"
#endif

#endif

