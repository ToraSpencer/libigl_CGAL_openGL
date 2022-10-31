#include "collapse_edge.h"
#include "circulation.h"
#include "edge_collapse_is_valid.h"
#include "decimate_trivial_callbacks.h"
#include <vector>


// 重载1.1——折叠一条边；
IGL_INLINE bool igl::collapse_edge(
  const int uEdgeIdx, 
  const Eigen::RowVectorXd & collapsedVer, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  int & e1, 
  int & e2, 
  int & f1, 
  int & f2)
{
  std::vector<int>  nbrTrisIdx_src, nbrVersIdx_src; 
  circulation(uEdgeIdx,  true, tris, edgeUeInfo, EF, EI,  nbrVersIdx_src, nbrTrisIdx_src); 
  std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des; 
  circulation(uEdgeIdx,  false, tris, edgeUeInfo, EF, EI,  nbrVersIdx_des, nbrTrisIdx_des); 
  return collapse_edge(uEdgeIdx, collapsedVer, nbrVersIdx_src, nbrTrisIdx_src, nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, EF, EI, e1, e2, f1, f2); 
}


// 重载1——折叠一条边；
IGL_INLINE bool igl::collapse_edge(
  const int uEdgeIdx, 
  const Eigen::RowVectorXd & collapsedVer, 
  /*const*/ std::vector<int> & nbrVersIdx_src, 
  const std::vector<int> & nbrTrisIdx_src, 
  /*const*/ std::vector<int> & nbrVersIdx_des, 
  const std::vector<int> & nbrTrisIdx_des, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  int & a_e1, 
  int & a_e2, 
  int & a_f1, 
  int & a_f2)
{
    /*
       Assign this to 0 rather than,  say,  -1 so that deleted elements will get draw as degenerate elements at vertex 0 
              (which should always exist and never get collapsed to anything else since it is the smallest index)
   */
  using namespace Eigen; 
  using namespace std; 
  const int eFlipFlag = uEdges(uEdgeIdx, 0) > uEdges(uEdgeIdx, 1); 

  // 定义无向边的起点和终点——索引小的为起点，索引大的为终点；
  const int srcIdx = eFlipFlag ? uEdges(uEdgeIdx, 1) : uEdges(uEdgeIdx, 0); 
  const int desIdx = eFlipFlag ? uEdges(uEdgeIdx, 0) : uEdges(uEdgeIdx, 1); 

  if(!edge_collapse_is_valid(nbrVersIdx_src, nbrVersIdx_des))
        return false; 

  // Important to grab neighbors of desIdx before monkeying with uEdges
  const std::vector<int>& nV2Fd = (!eFlipFlag ? nbrTrisIdx_src : nbrTrisIdx_des); 

  assert(srcIdx<desIdx && "srcIdx should be less than desIdx");  // The following implementation strongly relies on srcIdx<desIdx

  // move source and destination to placement
  vers.row(srcIdx) = collapsedVer; 
  vers.row(desIdx) = collapsedVer; 

  // lambda——replace edge and associate information with NULL
  const auto & kill_edge = [&uEdges, &EI, &EF](const int uEdgeIdx)
  {
    uEdges(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL; 
    uEdges(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL; 
    EF(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL; 
    EF(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL; 
    EI(uEdgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL; 
    EI(uEdgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL; 
  }; 

  // update edge info for each flap
  const int trisCount = tris.rows(); 
  for(int side = 0; side<2; side++)
  {
    const int f = EF(uEdgeIdx, side); 
    const int v = EI(uEdgeIdx, side); 
    const int sign = (eFlipFlag==0?1:-1)*(1-2*side); 
    // next edge emanating from desIdx
    const int e1 = edgeUeInfo(f+trisCount*((v+sign*1+3)%3)); 
    // prev edge pointing to srcIdx
    const int e2 = edgeUeInfo(f+trisCount*((v+sign*2+3)%3)); 
    assert(uEdges(e1, 0) == desIdx || uEdges(e1, 1) == desIdx); 
    assert(uEdges(e2, 0) == srcIdx || uEdges(e2, 1) == srcIdx); 
    // face adjacent to f on e1,  also incident on desIdx
    const bool flip1 = EF(e1, 1)==f; 
    const int f1 = flip1 ? EF(e1, 0) : EF(e1, 1); 
    assert(f1!=f); 
    assert(tris(f1, 0)==desIdx || tris(f1, 1)==desIdx || tris(f1, 2) == desIdx); 
    // across from which vertex of f1 does e1 appear?
    const int v1 = flip1 ? EI(e1, 0) : EI(e1, 1); 
    // Kill e1
    kill_edge(e1); 
    // Kill f
    tris(f, 0) = IGL_COLLAPSE_EDGE_NULL; 
    tris(f, 1) = IGL_COLLAPSE_EDGE_NULL; 
    tris(f, 2) = IGL_COLLAPSE_EDGE_NULL; 
    // map f1'srcIdx edge on e1 to e2
    assert(edgeUeInfo(f1+trisCount*v1) == e1); 
    edgeUeInfo(f1+trisCount*v1) = e2; 
    // side opposite f2,  the face adjacent to f on e2,  also incident on srcIdx
    const int opp2 = (EF(e2, 0)==f?0:1); 
    assert(EF(e2, opp2) == f); 
    EF(e2, opp2) = f1; 
    EI(e2, opp2) = v1; 
    // remap e2 from desIdx to srcIdx
    uEdges(e2, 0) = uEdges(e2, 0)==desIdx ? srcIdx : uEdges(e2, 0); 
    uEdges(e2, 1) = uEdges(e2, 1)==desIdx ? srcIdx : uEdges(e2, 1); 
    if(side==0)
    {
      a_e1 = e1; 
      a_f1 = f; 
    }else
    {
      a_e2 = e1; 
      a_f2 = f; 
    }
  }

  /*
       finally,  reindex faces and uEdges incident on desIdx. Do this last so asserts make sense.
       Could actually skip first and last,  since those are always the two collpased faces. 
       Nah,  this is handled by (tris(f, v) == desIdx)
       Don't attempt to use Nde, Nse here because edgeUeInfo has changed
  */
  {
    int p1 = -1; 
    for(auto f : nV2Fd)
    {
      for(int v = 0; v<3; v++)
      {
        if(tris(f, v) == desIdx)
        {
          const int e1 = edgeUeInfo(f+trisCount*((v+1)%3)); 
          const int flip1 = (EF(e1, 0)==f)?1:0; 
          assert( uEdges(e1, flip1) == desIdx || uEdges(e1, flip1) == srcIdx); 
          uEdges(e1, flip1) = srcIdx; 
          const int e2 = edgeUeInfo(f+trisCount*((v+2)%3)); 

          // Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
          if(e2 != p1)
          {
            const int flip2 = (EF(e2, 0)==f)?0:1; 
            assert( uEdges(e2, flip2) == desIdx || uEdges(e2, flip2) == srcIdx); 
            uEdges(e2, flip2) = srcIdx; 
          }

          tris(f, v) = srcIdx; 
          p1 = e1; 
          break; 
        }
      }
    }
  }


  // Finally,  "remove" this edge and its information
  kill_edge(uEdgeIdx); 
  return true; 
}


// 重载1.2
IGL_INLINE bool igl::collapse_edge(
  const int uEdgeIdx, 
  const Eigen::RowVectorXd & collapsedVer, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI)
{
  int e1, e2, f1, f2; 
  return collapse_edge(uEdgeIdx, collapsedVer, vers, tris, uEdges, edgeUeInfo, EF, EI, e1, e2, f1, f2); 
}


// 重载2.1
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & timeStamps, 
  Eigen::MatrixXd & collapsedVers)
{
  int uEdgeIdx, e1, e2, f1, f2; 
  decimate_pre_collapse_callback always_try; 
  decimate_post_collapse_callback never_care; 

  // 获取默认的pre_collapse(), post_collapse()函数子；
  decimate_trivial_callbacks(always_try, never_care); 
  
  return 
    collapse_edge(cost_and_placement, always_try, never_care, \
            vers, tris, uEdges, edgeUeInfo, EF, EI, workQueue, timeStamps, collapsedVers, uEdgeIdx, e1, e2, f1, f2); 
}


// 重载2.2
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_pre_collapse_callback       & pre_collapse, 
  const decimate_post_collapse_callback      & post_collapse, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & timeStamps, 
  Eigen::MatrixXd & collapsedVers)
{
  int uEdgeIdx, e1, e2, f1, f2; 
  return 
    collapse_edge(
      cost_and_placement, pre_collapse, post_collapse, 
      vers, tris, uEdges, edgeUeInfo, EF, EI, workQueue, timeStamps, collapsedVers, uEdgeIdx, e1, e2, f1, f2); 
}


// 重载2——循环调用折叠单调边的collapse_edge()重载1，直到满足终止条件；
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_pre_collapse_callback       & pre_collapse, 
  const decimate_post_collapse_callback      & post_collapse, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & timeStamps, 
  Eigen::MatrixXd & collapsedVers, 
  int & uEdgeIdx, 
  int & e1, 
  int & e2, 
  int & f1, 
  int & f2)
{
  using namespace Eigen; 
  using namespace igl; 

  // 1. 队列的循环
  std::tuple<double, int, int> edgeTuple; 
  while(true)
  {
    if(workQueue.empty())   // no uEdges to collapse
    {
      uEdgeIdx = -1; 
      return false; 
    }

    // 出队：
    edgeTuple = workQueue.top(); 
    if(std::get<0>(edgeTuple) == std::numeric_limits<double>::infinity())
    {
      uEdgeIdx = -1;      // min cost edge is infinite cost
      return false; 
    }
    workQueue.pop(); 
    uEdgeIdx = std::get<1>(edgeTuple); 

    // Check if matches timestamp
    if(std::get<2>(edgeTuple) == timeStamps(uEdgeIdx))
         break; 

    assert(std::get<2>(edgeTuple)  < timeStamps(uEdgeIdx) || timeStamps(uEdgeIdx) == -1);          // must be stale or dead.
  }
 
  // If we just need original face neighbors of edge,  could we gather that more directly than gathering face neighbors of each vertex?
  
  // 2. 计算当前边两端点1领域的顶点、三角片：
  std::vector<int> nbrTrisIdx_src, nbrVersIdx_src; 
  circulation(uEdgeIdx, true, tris, edgeUeInfo, EF, EI,  nbrVersIdx_src, nbrTrisIdx_src); 
  std::vector<int>  nbrTrisIdx_des, nbrVersIdx_des; 
  circulation(uEdgeIdx,  false, tris, edgeUeInfo, EF, EI,  nbrVersIdx_des, nbrTrisIdx_des); 


  // 3. 折叠队首的边：
  bool collapsed = true; 
  if(pre_collapse(vers, tris, uEdges, edgeUeInfo, EF, EI, workQueue, timeStamps, collapsedVers, uEdgeIdx)) // 默认情形下pre_collapse()什么都不做
        collapsed = collapse_edge(uEdgeIdx, collapsedVers.row(uEdgeIdx), nbrVersIdx_src, nbrTrisIdx_src, nbrVersIdx_des, nbrTrisIdx_des, vers, tris, uEdges, edgeUeInfo, EF, EI, e1, e2, f1, f2); 
  else
        collapsed = false;                       // Aborted by pre collapse callback

  post_collapse(vers, tris, uEdges, edgeUeInfo, EF, EI, workQueue, timeStamps, \
        collapsedVers, uEdgeIdx, e1, e2, f1, f2, collapsed);          // 默认情形下post_collapse()什么都不做

   // 4. 
  if(collapsed)
  {
    // 4.1 Erase the two,  other collapsed uEdges by marking their timestamps as -1
    timeStamps(e1) = -1; 
    timeStamps(e2) = -1; 

    // TODO: visits uEdges multiple times,  ~150% more updates than should
    /*
         update local neighbors
         loop over original face neighbors
    
         Can't use previous computed Nse and Nde because those refer to edgeUeInfo
         before it was changed...
    */

    // 4.2
    std::vector<int> nbrTrisIdx; 
    nbrTrisIdx.reserve( nbrTrisIdx_src.size() + nbrTrisIdx_des.size() );                           // preallocate memory
    nbrTrisIdx.insert( nbrTrisIdx.end(),  nbrTrisIdx_src.begin(),  nbrTrisIdx_src.end() ); 
    nbrTrisIdx.insert( nbrTrisIdx.end(),  nbrTrisIdx_des.begin(),  nbrTrisIdx_des.end() ); 
    std::sort( nbrTrisIdx.begin(),  nbrTrisIdx.end() );                               // https://stackoverflow.com/a/1041939/148668
    nbrTrisIdx.erase( std::unique( nbrTrisIdx.begin(),  nbrTrisIdx.end() ),  nbrTrisIdx.end() ); 

    // 4.3 Collect all uEdges that must be updated
    std::vector<int> Ne; 
    Ne.reserve(3 * nbrTrisIdx.size()); 
    for(auto & triIdx : nbrTrisIdx)
    {
      if(tris(triIdx, 0) != IGL_COLLAPSE_EDGE_NULL ||
          tris(triIdx, 1) != IGL_COLLAPSE_EDGE_NULL ||
          tris(triIdx, 2) != IGL_COLLAPSE_EDGE_NULL)
      {
        for(int i = 0;  i < 3;  i++)
        {         
          const int ueIdx = edgeUeInfo(i * tris.rows() + triIdx);          
          Ne.push_back(ueIdx); 
        }
      }
    }

    // Only process edge once
    std::sort( Ne.begin(),  Ne.end() ); 
    Ne.erase( std::unique( Ne.begin(),  Ne.end() ),  Ne.end() );             // 去除重复元素：
    for(auto & ueIdx : Ne)
    {
       // 计算边折叠的cost值，及折叠后的顶点坐标：
       double cost; 
       RowVectorXd place; 
       cost_and_placement(ueIdx, vers, tris, uEdges, edgeUeInfo, EF, EI, cost, place); 

       // Increment timestamp
       timeStamps(ueIdx)++; 

       // Replace in queue
       workQueue.emplace(cost, ueIdx, timeStamps(ueIdx)); 
       collapsedVers.row(ueIdx) = place; 
    }
  }
  else
  {
    // reinsert with infinite weight (the provided cost function must **not** have given this un-collapsable edge inf cost already)
    timeStamps(uEdgeIdx)++;          // Increment timestamp

    // Replace in queue
    workQueue.emplace(std::numeric_limits<double>::infinity(), uEdgeIdx, timeStamps(uEdgeIdx)); 
  }


  return collapsed; 
}
