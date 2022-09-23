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
IGL_INLINE void igl::bfs(const std::vector<std::vector<AType> > & adjList,           // �ڽӱ�
      const size_t startIdx,
      std::vector<DType> & disCoveredIdx,                // ������ȱ���֮�����б����ʹ��Ķ����������
      std::vector<PType> & bfsTreeVec)          // ��ʾ���������������������������������Ϊi�Ľ�㣬�丸�ڵ�����ΪbfsTreeVec[i]
{
      int versCount = startIdx+1;
      for(const auto & vec : adjList)
          for(const auto & index : vec) 
              versCount = std::max(versCount, index+1);
      
      std::vector<bool> visited(versCount, false);              // ��Ƕ����Ƿ��ѱ����ʣ�
      bfsTreeVec.resize(versCount, -1);
      std::queue<std::pair<int, int>> workingQueue;
      workingQueue.push({startIdx, -1});

      while(!workingQueue.empty())
      {
        const int index1 = workingQueue.front().first;                  // ��ǰ���㣻
        const int index2 = workingQueue.front().second;             // ��ǰ����1�������ѷ��ʹ��Ķ���
        workingQueue.pop();     // ��ǰ������ӣ�
        if(visited[index1])     
            continue;

        disCoveredIdx.push_back(index1);
        bfsTreeVec[index1] = index2;
        visited[index1] = true;                                  // ��ǰ������Ϊ�ѷ���״̬��
        for(const auto & adjIdx1 : adjList[index1]) 
            workingQueue.push({adjIdx1, index1});   // ��ǰ����1������û�з��ʹ��Ķ�����ӣ�
      }
}


// imp3:
template <
  typename AType,
  typename DType,
  typename PType>
IGL_INLINE void igl::bfs(
  const Eigen::SparseCompressedBase<AType> & A,         // �ڽӾ���
  const size_t startIdx,
  std::vector<DType> & disCoveredIdx,           // ������ȱ���֮�����б����ʹ��Ķ����������
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

