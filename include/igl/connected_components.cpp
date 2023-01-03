#include "connected_components.h"
#include <queue>

template < typename Atype, typename DerivedC, typename DerivedK>
IGL_INLINE int igl::connected_components(
  const Eigen::SparseMatrix<Atype> & adjSM,
  Eigen::PlainObjectBase<DerivedC> & connectedLabels,
  Eigen::PlainObjectBase<DerivedK> & connectedCount)
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

  typedef typename Eigen::SparseMatrix<Atype>::Index Index;
  const auto versCount = adjSM.rows();
  assert(adjSM.cols() == adjSM.rows() && "adjSM should be square");

  // 1. 初始化：
  connectedLabels.setConstant(versCount,1,versCount);           // versCount  means not yet visited
  connectedCount.setZero(versCount,1);

  // 2. 遍历所有顶点，搜索与其联通的顶点：
  typename DerivedC::Scalar currentLabel = 0;
  for(Eigen::Index i = 0; i < versCount; i++)
  {
    // 2.1 若当前顶点已被访问过，continue:
    if(connectedLabels(i) < versCount) 
        continue;

    // 2.2 若当前顶点未被访问过，执行BFS
    std::queue<Index> workingQueue;
    workingQueue.push(i);
    while(!workingQueue.empty())
    {
      const Index curVerIdx = workingQueue.front();
      workingQueue.pop();

      if(connectedLabels(curVerIdx) < versCount) 
          continue;

      //    当前顶点label赋值，标记为被访问；
      connectedLabels(curVerIdx) = currentLabel;
      connectedCount(currentLabel)++;

      //     在邻接矩阵中搜索当前顶点邻接的顶点：
      for(typename Eigen::SparseMatrix<Atype>::InnerIterator iter (adjSM, curVerIdx); iter; ++iter)
      {
        const Index connectVerIdx = iter.index();

        if(connectedLabels(connectVerIdx) < versCount) 
            continue;

        workingQueue.push(connectVerIdx);
      }
    }


    // 2.3 上一个标签的顶点收集完毕，下一个循环收集下一个标签的顶点：
    currentLabel++;
  }

  // 3. shrink_to_fit()
  connectedCount.conservativeResize(currentLabel,1);

  return currentLabel;
}

#ifdef IGL_STATIC_LIBRARY
template int igl::connected_components<int, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<int, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
