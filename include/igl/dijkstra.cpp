#include "dijkstra.h"


// dijkstra()����1
template <typename IndexType, typename DerivedD, typename DerivedP>
IGL_INLINE int igl::dijkstra(
  const IndexType &startIdx,
  const std::set<IndexType> &targets,
  const std::vector<std::vector<IndexType> >& adjList,
  const std::vector<double>& weights,
  Eigen::PlainObjectBase<DerivedD> &minDis,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
  int versCount = adjList.size();
  minDis.setConstant(versCount, 1, std::numeric_limits<typename DerivedD::Scalar>::max());
  minDis[startIdx] = 0;
  previous.setConstant(versCount, 1, -1);
  std::set<std::pair<typename DerivedD::Scalar, IndexType> > vertex_queue;
  vertex_queue.insert(std::make_pair(minDis[startIdx], startIdx));

  while (!vertex_queue.empty())
  {
    typename DerivedD::Scalar dist = vertex_queue.begin()->first;
    IndexType index0 = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());

    if (targets.find(index0)!= targets.end())
      return index0;

    // Visit each edge exiting index0
    const std::vector<int> &adjVersIdx0 = adjList[index0];
    for (std::vector<int>::const_iterator adjIter = adjVersIdx0.begin();
         adjIter != adjVersIdx0.end();
         adjIter++)
    {
      IndexType adjIdx0 = *adjIter;
      typename DerivedD::Scalar pathLen = dist + weights[index0];
      if (pathLen < minDis[adjIdx0])
      {
        vertex_queue.erase(std::make_pair(minDis[adjIdx0], adjIdx0));

        minDis[adjIdx0] = pathLen;
        previous[adjIdx0] = index0;
        vertex_queue.insert(std::make_pair(minDis[adjIdx0], adjIdx0));
      }
    }
  }
  //we should never get here
  return -1;
}


// dijkstra()����1.1
template <typename IndexType, typename DerivedD, typename DerivedP>
IGL_INLINE int igl::dijkstra(
  const IndexType &startIdx,
  const std::set<IndexType> &targets,
  const std::vector<std::vector<IndexType> >& adjList,
  Eigen::PlainObjectBase<DerivedD> &minDis,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
  std::vector<double> weights(adjList.size(), 1.0);
  return dijkstra(startIdx, targets, adjList, weights, minDis, previous);
}


// dijkstra()����2����·��������дΪstd::vector��ʾ��
template <typename IndexType, typename DerivedP>
IGL_INLINE void igl::dijkstra(
  const IndexType &vertex,
  const Eigen::MatrixBase<DerivedP> &previous,
  std::vector<IndexType> &path)
{
  IndexType startIdx = vertex;
  path.clear();
  for ( ; startIdx != -1; startIdx = previous[startIdx])
    path.push_back(startIdx);
}


// dijkstra()����3
template <typename IndexType, typename DerivedV,
typename DerivedD, typename DerivedP>
IGL_INLINE int igl::dijkstra(
  const Eigen::MatrixBase<DerivedV> &vers,
  const std::vector<std::vector<IndexType> >& adjList,
  const IndexType &startIdx,
  const std::set<IndexType> &targets,
  Eigen::PlainObjectBase<DerivedD> &minDis,
  Eigen::PlainObjectBase<DerivedP> &previous)
{
    /*
    int igl::dijkstra(                                                                                                 ����·����̵��յ��������������յ�ʱ������-1
                const Eigen::MatrixBase<DerivedV> &vers,                                        ͼ�ĵ㼯
                const std::vector<std::vector<IndexType> >& adjList,                        �ڽӱ�
                const IndexType &startIdx,                                                                   �������
                const std::set<IndexType> &targets,                                                    �յ��������ɴ�������
                Eigen::PlainObjectBase<DerivedD> &minDis,        ��targetsΪ�գ���Ϊÿ�����㵽������·���ĳ�������������Ϊ�գ�����С������δ�������Ķ��㳤��Ϊinf   
                Eigen::PlainObjectBase<DerivedP> &previous      ��targetsΪ�գ���previousΪ��С����������������Ϊ�գ���Ϊ����������յ�ʱ����С��������
                  )
        
    */

  int versCount = adjList.size();

  minDis.setConstant(versCount, 1, std::numeric_limits<typename DerivedD::Scalar>::infinity());    // ��ʼ��С������Ϊinf
  minDis[startIdx] = 0;
  previous.setConstant(versCount, 1, -1);
  std::set<std::pair<typename DerivedD::Scalar, IndexType> > vertex_queue;      // working queue; (minDis, index)
  vertex_queue.insert(std::make_pair(minDis[startIdx], startIdx));

  // working queue�ڵ�ѭ����
  while (!vertex_queue.empty())
  {
    typename DerivedD::Scalar dist = vertex_queue.begin()->first;           // ����Ԫ�ض���index0��startIdx��С���룻
    IndexType index0 = vertex_queue.begin()->second;

    // 1. popHead();
    vertex_queue.erase(vertex_queue.begin());

    if (targets.find(index0)!= targets.end())
      return index0;

    // 2. ����index0��1���������ж��㣻����������index0���������бߣ�
    const std::vector<int> &adjVersIdx0 = adjList[index0];
    for (std::vector<int>::const_iterator adjIter = adjVersIdx0.begin(); adjIter != adjVersIdx0.end(); adjIter++)   
    {
      IndexType adjIdx0 = *adjIter;

      // 2.1 ���㾭����ǰ��(index0, adjIdx0)֮���ۼƵ�·�����ȣ�
      typename DerivedD::Scalar pathLen = dist + (vers.row(index0) - vers.row(adjIdx0)).norm();

      // 2.2 ��ʹ�õ�ǰ�ߣ���·�����ȱ�֮ǰ��С����ʹ�õ�ǰ�ߣ�
      if (pathLen < minDis[adjIdx0])
      {
        vertex_queue.erase(std::make_pair(minDis[adjIdx0], adjIdx0));

        minDis[adjIdx0] = pathLen;
        previous[adjIdx0] = index0;         
        vertex_queue.insert(std::make_pair(minDis[adjIdx0], adjIdx0));        // queue pushHead()
      }
    }
  }

  return -1;
}



// ģ���ػ���
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::dijkstra<int, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, std::set<int, std::less<int>, std::allocator<int> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template int igl::dijkstra<int, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, std::set<int, std::less<int>, std::allocator<int> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::dijkstra<int, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<int, std::allocator<int> >&);
template int igl::dijkstra<int, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, int const&, std::set<int, std::less<int>, std::allocator<int> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
