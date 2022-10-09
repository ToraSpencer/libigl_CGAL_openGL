#ifndef IGL_UNIQUE_EDGE_MAP_H
#define IGL_UNIQUE_EDGE_MAP_H
#include "igl_inline.h"
#include <Eigen/Dense>
#include <vector>


namespace igl
{
    // жиди1.1
    /*
       Construct relationships between facet "half"-(or rather "viewed")-edges edges
       to unique edges of the mesh seen as a graph.
  
       Inputs:
         tris           by 3  list of simplices

       Outputs:
         edges                    #tris*3 by 2 list of all directed edges, such that edges.row(f+#tris*c) is the edge opposite tris(f,c)
         uEdges                  #uEdges by 2 list of unique undirected edges
         edgeUeInfo           #tris*3 list of indices into uEdges, mapping each directed edge to unique
                                                undirected edge so that uEdges(edgeUeInfo(f+#tris*c)) is the unique edge
                                                corresponding to edges.row(f+#tris*c)
         uE2E         #uEdges list of lists of indices into edges of coexisting edges, so that
                                 edges.row(uE2E[i][j]) corresponds to uEdges.row(i) for all j in 0..uE2E[i].size()-1.
    */
  template <
    typename DerivedF,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2EType>
  IGL_INLINE void unique_edge_map(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DerivedE> & edges,
    Eigen::PlainObjectBase<DeriveduE> & uEdges,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo,
    std::vector<std::vector<uE2EType> > & uE2E);


  // жиди1.
  template <
    typename DerivedF,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP>
  IGL_INLINE void unique_edge_map(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DerivedE> & edges,
    Eigen::PlainObjectBase<DeriveduE> & uEdges,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo);

  // жиди1.2
  /*
   Outputs:
     uEC  #uEC+1 list of cumulative counts of directed edges sharing each
       unique edge so the uEC(i+1)-uEC(i) is the number of directed edges
       sharing the ith unique edge.

     uEE  #edges list of indices into edges, so that the consecutive segment of
       indices uEE.segment(uEC(i),uEC(i+1)-uEC(i)) lists all directed edges
       sharing the ith unique edge.
  */
  template <
    typename DerivedF,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename DeriveduEC,
    typename DeriveduEE>
  IGL_INLINE void unique_edge_map(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::PlainObjectBase<DerivedE> & edges,
    Eigen::PlainObjectBase<DeriveduE> & uEdges,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo,
    Eigen::PlainObjectBase<DeriveduEC> & uEC,
    Eigen::PlainObjectBase<DeriveduEE> & uEE);

}
#ifndef IGL_STATIC_LIBRARY
#  include "unique_edge_map.cpp"
#endif

#endif
