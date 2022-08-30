#ifndef IGL_BFS_H
#define IGL_BFS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <Eigen/Sparse>


// 图的广度优先搜索(BFS: breadth first search)
namespace igl
{
    /*
         Traverse a **directed** graph represented by an adjacency list using
       breadth first search
  
       Inputs:
         A  #V list of adjacency lists  or #V by #V adjacency matrix
         s  starting node (index into A)

       Outputs:
         disCoveredIdx  #V list of indices into rows of A in the order in which graph nodes
           are discovered.
         bfsTreeVec  #V list of indices into rows of A of predecessor in resulting
           spanning tree {-1 indicates root/not discovered), order corresponds to
           V **not** disCoveredIdx.
    */
  template <
    typename AType,
    typename DerivedD,
    typename DerivedP>
  IGL_INLINE void bfs(
    const AType & A,
    const size_t s,
    Eigen::PlainObjectBase<DerivedD> & disCoveredIdx,
    Eigen::PlainObjectBase<DerivedP> & bfsTreeVec);


  template <
    typename AType,
    typename DType,
    typename PType>
  IGL_INLINE void bfs(
    const std::vector<std::vector<AType> > & A,
    const size_t s,
    std::vector<DType> & disCoveredIdx,
    std::vector<PType> & bfsTreeVec);


  template <
    typename AType,
    typename DType,
    typename PType>
  IGL_INLINE void bfs(
    const Eigen::SparseCompressedBase<AType> & A,
    const size_t s,
    std::vector<DType> & disCoveredIdx,
    std::vector<PType> & bfsTreeVec);
}
#ifndef IGL_STATIC_LIBRARY
#  include "bfs.cpp"
#endif
#endif

