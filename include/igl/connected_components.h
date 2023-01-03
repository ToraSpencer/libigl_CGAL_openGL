#ifndef IGL_CONNECTED_COMPONENTS_H
#define IGL_CONNECTED_COMPONENTS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
    /*
       Determine the connected components of a graph described by the input adjacency matrix (similar to MATLAB's graphconncomp).
  
       Inputs:
          adjSM  #adjSM by #adjSM adjacency matrix (treated as describing an undirected graph)
                                     网格的无向边邻接矩阵；

       Outputs:
          connectedLabels  #adjSM list of component indices into [0,#connectedCount-1]
                                        同一连通区域内的顶点标签，依次为0, 1, 2,...
                                        索引为i的顶点的标签为connectedLabels(i);

          connectedCount  #connectedCount list of sizes of each component
                                        每个标签对应的单连通区域内包含的顶点数
                                        标签为i的单连通区域包含的顶点数为connectedCount(i)
       
       Returns number of connected components
    */
  template < typename Atype, typename DerivedC, typename DerivedK>
  IGL_INLINE int connected_components(
    const Eigen::SparseMatrix<Atype> & adjSM,
    Eigen::PlainObjectBase<DerivedC> & connectedLabels,
    Eigen::PlainObjectBase<DerivedK> & connectedCount);
}

#ifndef IGL_STATIC_LIBRARY
#  include "connected_components.cpp"
#endif

#endif
