#ifndef IGL_SHARP_EDGES_H
#define IGL_SHARP_EDGES_H

#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // SHARP_EDGES Given a mesh, compute sharp edges.
  //
  // Inputs:
  //   vers  #vers by 3 list of vertex positions
  //   tris  #tris by 3 list of triangle mesh indices into vers
  //   angle  dihedral angle considered to sharp (e.g., igl::PI * 0.11)
  // Outputs:
  //   SE  #SE by 2 list of edge indices into vers
  //   uE  #uE by 2 list of unique undirected edges
  //   EMAP #tris*3 list of indices into uE, mapping each directed edge to unique
  //     undirected edge so that uE(EMAP(f+#tris*c)) is the unique edge
  //     corresponding to E.row(f+#tris*c)
  //   uE2E  #uE list of lists of indices into E of coexisting edges, so that
  //     E.row(uE2E[i][j]) corresponds to uE.row(i) for all j in
  //     0..uE2E[i].size()-1.
  //   sharp  #SE list of indices into uE revealing sharp undirected edges
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedSE,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2Etype,
    typename sharptype>
  IGL_INLINE void sharp_edges(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const typename DerivedV::Scalar angle,
    Eigen::PlainObjectBase<DerivedSE> & SE,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DeriveduE> & uE,
    Eigen::PlainObjectBase<DerivedEMAP> & EMAP,
    std::vector<std::vector<uE2Etype> > & uE2E,
    std::vector< sharptype > & sharp);
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedSE>
  IGL_INLINE void sharp_edges(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const typename DerivedV::Scalar angle,
    Eigen::PlainObjectBase<DerivedSE> & SE
    );
}

#ifndef IGL_STATIC_LIBRARY
#  include "sharp_edges.cpp"
#endif

#endif 
