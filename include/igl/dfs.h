#ifndef IGL_DFS_H
#define IGL_DFS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

// 图的深度优先搜索(DFS: depth first search)
namespace igl
{
    /*
       Traverse a **directed** graph represented by an adjacency list using
       depth first search
  
       Inputs:
         adjList  #vers list of adjacency lists
         startIdx  starting node (index into adjList)

       Outputs:
         discoveredIdx  #vers list of indices into rows of adjList in the order in which graph nodes are discovered.

         dfsTreeVec  #vers list of indices into rows of adjList of predecessor in resulting
           spanning tree {-1 indicates root/not discovered), order corresponds to
           vers **not** discoveredIdx.

         closedIdx  #vers list of indices into rows of adjList in order that nodes are "closed"
           (all descendants have been discovered)
  */
  template <typename AType,
    typename DerivedD,
    typename DerivedP,
    typename DerivedC>
  IGL_INLINE void dfs(
    const std::vector<std::vector<AType> > & adjList,
    const size_t startIdx,
    Eigen::PlainObjectBase<DerivedD> & discoveredIdx,
    Eigen::PlainObjectBase<DerivedP> & dfsTreeVec,
    Eigen::PlainObjectBase<DerivedC> & closedIdx);


  template <
    typename AType,
    typename DType,
    typename PType,
    typename CType>
  IGL_INLINE void dfs(
    const std::vector<std::vector<AType> > & adjList,
    const size_t startIdx,
    std::vector<DType> & discoveredIdx,
    std::vector<PType> & dfsTreeVec,
    std::vector<CType> & closedIdx);

}
#ifndef IGL_STATIC_LIBRARY
#  include "dfs.cpp"
#endif
#endif
