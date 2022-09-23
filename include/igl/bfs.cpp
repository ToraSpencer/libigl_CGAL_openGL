#include "bfs.h"
#include "list_to_matrix.h"
#include <vector>
#include <queue>

// imp1
template <typename AType,
  typename DerivedD,
  typename DerivedP>
IGL_INLINE void igl::bfs(const AType & A,  const size_t startIdx,
      Eigen::PlainObjectBase<DerivedD> & disCoveredIdx,
      Eigen::PlainObjectBase<DerivedP> & bfsTreeVec)
{
  std::vector<typename DerivedD::Scalar> vD;
  std::vector<typename DerivedP::Scalar> vP;
  bfs(A,startIdx,vD,vP);
  list_to_matrix(vD,disCoveredIdx);
  list_to_matrix(vP, bfsTreeVec);
}


// imp2:
template <typename AType,  typename DType,  typename PType>
IGL_INLINE void igl::bfs(const std::vector<std::vector<AType> > & adjList,           // 邻接表
      const size_t startIdx,
      std::vector<DType> & disCoveredIdx,                // 广度优先遍历之后所有被访问过的顶点的索引；
      std::vector<PType> & bfsTreeVec)          // 表示广度优先搜索树的索引向量；树中索引为i的结点，其父节点索引为bfsTreeVec[i]
{
      int versCount = startIdx+1;
      for(const auto & vec : adjList)
          for(const auto & index : vec) 
              versCount = std::max(versCount, index+1);
      
      std::vector<bool> visited(versCount, false);              // 标记顶点是否已被访问；
      bfsTreeVec.resize(versCount, -1);
      std::queue<std::pair<int, int>> workingQueue;
      workingQueue.push({startIdx, -1});

      while(!workingQueue.empty())
      {
        const int index1 = workingQueue.front().first;                  // 当前顶点；
        const int index2 = workingQueue.front().second;             // 当前顶点1领域内已访问过的顶点
        workingQueue.pop();     // 当前顶点出队；
        if(visited[index1])     
            continue;

        disCoveredIdx.push_back(index1);
        bfsTreeVec[index1] = index2;
        visited[index1] = true;                                  // 当前顶点标记为已访问状态；
        for(const auto & adjIdx1 : adjList[index1]) 
            workingQueue.push({adjIdx1, index1});   // 当前顶点1领域内没有访问过的顶点入队；
      }
}


// imp3:
template <
  typename AType,
  typename DType,
  typename PType>
IGL_INLINE void igl::bfs(
  const Eigen::SparseCompressedBase<AType> & A,         // 邻接矩阵
  const size_t startIdx,
  std::vector<DType> & disCoveredIdx,           // 广度优先遍历之后所有被访问过的顶点的索引；
  std::vector<PType> & bfsTreeVec)
{
  // number of nodes
  int versCount = A.rows();
  assert(A.rows() == A.cols());
  std::vector<bool> visited(versCount,false);
  bfsTreeVec.resize(versCount,-1);
  std::queue<std::pair<int,int> > workingQueue;
  workingQueue.push({startIdx,-1});

  while(!workingQueue.empty())
  {
    const int index1 = workingQueue.front().first;
    const int index2 = workingQueue.front().second;
    workingQueue.pop();

    if(visited[index1])
      continue;

    disCoveredIdx.push_back(index1);
    bfsTreeVec[index1] = index2;
    visited[index1] = true;
    for(typename AType::InnerIterator it (A,index1); it; ++it)
      if(it.value()) 
          workingQueue.push({it.index(),index1});
  }
}


#ifdef IGL_STATIC_LIBRARY
template void igl::bfs<std::vector<std::vector<int, std::allocator<int>>, std::allocator<std::vector<int, std::allocator<int>>>>,\
    Eigen::Matrix<int, -1, 1, 0, -1, 1>, \
    Eigen::Matrix<int, -1, 1, 0, -1, 1>>\
    (   std::vector<std::vector<int, std::allocator<int>>, std::allocator<std::vector<int, std::allocator<int> > > > const&, \
        const size_t, \
        Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, \
        Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
 
#endif

