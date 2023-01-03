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
                                 �����������ڽӾ���

   Outputs:
      connectedLabels  #adjSM list of component indices into [0,#connectedCount-1]
                                    ͬһ��ͨ�����ڵĶ����ǩ������Ϊ0, 1, 2,...
                                    ����Ϊi�Ķ���ı�ǩΪconnectedLabels(i);

      connectedCount  #connectedCount list of sizes of each component
                                    ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ����Ķ�����
                                    ��ǩΪi�ĵ���ͨ��������Ķ�����ΪconnectedCount(i)

   Returns number of connected components
*/

  typedef typename Eigen::SparseMatrix<Atype>::Index Index;
  const auto versCount = adjSM.rows();
  assert(adjSM.cols() == adjSM.rows() && "adjSM should be square");

  // 1. ��ʼ����
  connectedLabels.setConstant(versCount,1,versCount);           // versCount  means not yet visited
  connectedCount.setZero(versCount,1);

  // 2. �������ж��㣬����������ͨ�Ķ��㣺
  typename DerivedC::Scalar currentLabel = 0;
  for(Eigen::Index i = 0; i < versCount; i++)
  {
    // 2.1 ����ǰ�����ѱ����ʹ���continue:
    if(connectedLabels(i) < versCount) 
        continue;

    // 2.2 ����ǰ����δ�����ʹ���ִ��BFS
    std::queue<Index> workingQueue;
    workingQueue.push(i);
    while(!workingQueue.empty())
    {
      const Index curVerIdx = workingQueue.front();
      workingQueue.pop();

      if(connectedLabels(curVerIdx) < versCount) 
          continue;

      //    ��ǰ����label��ֵ�����Ϊ�����ʣ�
      connectedLabels(curVerIdx) = currentLabel;
      connectedCount(currentLabel)++;

      //     ���ڽӾ�����������ǰ�����ڽӵĶ��㣺
      for(typename Eigen::SparseMatrix<Atype>::InnerIterator iter (adjSM, curVerIdx); iter; ++iter)
      {
        const Index connectVerIdx = iter.index();

        if(connectedLabels(connectVerIdx) < versCount) 
            continue;

        workingQueue.push(connectVerIdx);
      }
    }


    // 2.3 ��һ����ǩ�Ķ����ռ���ϣ���һ��ѭ���ռ���һ����ǩ�Ķ��㣺
    currentLabel++;
  }

  // 3. shrink_to_fit()
  connectedCount.conservativeResize(currentLabel,1);

  return currentLabel;
}

#ifdef IGL_STATIC_LIBRARY
template int igl::connected_components<int, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<int, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
