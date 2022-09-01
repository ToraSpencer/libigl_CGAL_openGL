#include "dijkstra.h"


// dijkstra()重载1
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


// dijkstra()重载1.1
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


// dijkstra()重载2——路径向量改写为std::vector表示；
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


// dijkstra()重载3
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
    int igl::dijkstra(                                                                                                 返回路径最短的终点索引；不输入终点时，返回-1
                const Eigen::MatrixBase<DerivedV> &vers,                                        图的点集
                const std::vector<std::vector<IndexType> >& adjList,                        邻接表
                const IndexType &startIdx,                                                                   起点索引
                const std::set<IndexType> &targets,                                                    终点索引，可传入多个；
                Eigen::PlainObjectBase<DerivedD> &minDis,        若targets为空，则为每个顶点到起点最短路径的长度向量；若不为空，则最小生成树未生长到的顶点长度为inf   
                Eigen::PlainObjectBase<DerivedP> &previous      若targets为空，则previous为最小生成树向量；若不为空，则为生长到最近终点时的最小生成树；
                  )
        
    */

  int versCount = adjList.size();

  minDis.setConstant(versCount, 1, std::numeric_limits<typename DerivedD::Scalar>::infinity());    // 初始最小距离设为inf
  minDis[startIdx] = 0;
  previous.setConstant(versCount, 1, -1);
  std::set<std::pair<typename DerivedD::Scalar, IndexType> > vertex_queue;      // working queue; (minDis, index)
  vertex_queue.insert(std::make_pair(minDis[startIdx], startIdx));

  // working queue内的循环：
  while (!vertex_queue.empty())
  {
    typename DerivedD::Scalar dist = vertex_queue.begin()->first;           // 队首元素顶点index0到startIdx最小距离；
    IndexType index0 = vertex_queue.begin()->second;

    // 1. popHead();
    vertex_queue.erase(vertex_queue.begin());

    if (targets.find(index0)!= targets.end())
      return index0;

    // 2. 遍历index0的1领域内所有顶点；即遍历顶点index0关联的所有边；
    const std::vector<int> &adjVersIdx0 = adjList[index0];
    for (std::vector<int>::const_iterator adjIter = adjVersIdx0.begin(); adjIter != adjVersIdx0.end(); adjIter++)   
    {
      IndexType adjIdx0 = *adjIter;

      // 2.1 计算经过当前边(index0, adjIdx0)之后累计的路径长度；
      typename DerivedD::Scalar pathLen = dist + (vers.row(index0) - vers.row(adjIdx0)).norm();

      // 2.2 若使用当前边，总路径长度比之前的小，则使用当前边；
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



// 模板特化：
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::dijkstra<int, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, std::set<int, std::less<int>, std::allocator<int> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template int igl::dijkstra<int, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, std::set<int, std::less<int>, std::allocator<int> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::dijkstra<int, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<int, std::allocator<int> >&);
template int igl::dijkstra<int, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, int const&, std::set<int, std::less<int>, std::allocator<int> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
