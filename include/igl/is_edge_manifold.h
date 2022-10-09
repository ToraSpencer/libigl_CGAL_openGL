#ifndef IGL_IS_EDGE_MANIFOLD_H
#define IGL_IS_EDGE_MANIFOLD_H
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl 
{
  // 检测网格每条边是否都是流形边
  /*
        每条流形边应该关联
      every edge is incident one one face (boundary) or two oppositely oriented faces 
  
       Inputs:
         tris  #tris by 3 list of triangle indices

       Returns true iff all edges are manifold
  
       See also: is_vertex_manifold
    */
  template <typename DerivedF>
  IGL_INLINE bool is_edge_manifold(
    const Eigen::MatrixBase<DerivedF>& tris);


  /*
       Inputs:
         tris                #tris by 3 list of triangle indices

       Outputs:
         BF                  #tris by 3 list of flags revealing if edge opposite corresponding vertex  is non-manifold.
         uEdges         #uEdges by 2 list of unique edges
         EMAP           3*#tris list of indices of opposite edges in "uEdges"
         BE                 #uEdges list of flages whether edge is non-manifold
  */
  template <
    typename DerivedF, 
    typename DerivedBF,
    typename DerivedE,
    typename DerivedEMAP,
    typename DerivedBE>
  IGL_INLINE bool is_edge_manifold(
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedBF>& BF,
    Eigen::PlainObjectBase<DerivedE>& uEdges,
    Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
    Eigen::PlainObjectBase<DerivedBE>& BE);


  // Inputs:
  //   tris  #tris by 3 list of triangle indices
  //   ne  number of edges (#uEdges)
  //   EMAP  3*#tris list of indices of opposite edges in "uEdges"
  // Outputs:
  //   BF  #tris by 3 list of flags revealing if edge opposite corresponding vertex
  //     is non-manifold.
  //   BE  ne list of flages whether edge is non-manifold
  template <
    typename DerivedF,
    typename DerivedEMAP,
    typename DerivedBF,
    typename DerivedBE>
  IGL_INLINE bool is_edge_manifold(
    const Eigen::MatrixBase<DerivedF>& tris,
    const typename DerivedF::Index ne,
    const Eigen::MatrixBase<DerivedEMAP>& EMAP,
    Eigen::PlainObjectBase<DerivedBF>& BF,
    Eigen::PlainObjectBase<DerivedBE>& BE);
}

#ifndef IGL_STATIC_LIBRARY
#  include "is_edge_manifold.cpp"
#endif

#endif
