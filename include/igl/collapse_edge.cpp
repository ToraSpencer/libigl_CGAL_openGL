#include "collapse_edge.h"
#include "circulation.h"
#include "edge_collapse_is_valid.h"
#include "decimate_trivial_callbacks.h"
#include <vector>


// 重载1.1
IGL_INLINE bool igl::collapse_edge(
  const int edgeIdx, 
  const Eigen::RowVectorXd & p, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  int & e1, 
  int & e2, 
  int & f1, 
  int & f2)
{
  std::vector<int> /*Nse, */Nsf, Nsv;
  circulation(edgeIdx,  true, tris, EMAP, EF, EI, /*Nse, */Nsv, Nsf);
  std::vector<int> /*Nde, */Ndf, Ndv;
  circulation(edgeIdx,  false, tris, EMAP, EF, EI, /*Nde, */Ndv, Ndf);
  return collapse_edge(edgeIdx, p, Nsv, Nsf, Ndv, Ndf, vers, tris, edges, EMAP, EF, EI, e1, e2, f1, f2);
}


// 重载1
IGL_INLINE bool igl::collapse_edge(
  const int edgeIdx, 
  const Eigen::RowVectorXd & p, 
  /*const*/ std::vector<int> & Nsv, 
  const std::vector<int> & Nsf, 
  /*const*/ std::vector<int> & Ndv, 
  const std::vector<int> & Ndf, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  int & a_e1, 
  int & a_e2, 
  int & a_f1, 
  int & a_f2)
{
  // Assign this to 0 rather than,  say,  -1 so that deleted elements will get
  // draw as degenerate elements at vertex 0 (which should always exist and
  // never get collapsed to anything else since it is the smallest index)
  using namespace Eigen;
  using namespace std;
  const int eflip = edges(edgeIdx, 0)>edges(edgeIdx, 1);
  // source and destination
  const int s = eflip?edges(edgeIdx, 1):edges(edgeIdx, 0);
  const int d = eflip?edges(edgeIdx, 0):edges(edgeIdx, 1);

  if(!edge_collapse_is_valid(Nsv, Ndv))
        return false;

  // OVERLOAD: caller may have just computed this
  //
  // Important to grab neighbors of d before monkeying with edges
  const std::vector<int> & nV2Fd = (!eflip ? Nsf : Ndf);

  // The following implementation strongly relies on s<d
  assert(s<d && "s should be less than d");
  // move source and destination to placement
  vers.row(s) = p;
  vers.row(d) = p;

  // Helper function to replace edge and associate information with NULL
  const auto & kill_edge = [&edges, &EI, &EF](const int edgeIdx)
  {
    edges(edgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
    edges(edgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
    EF(edgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
    EF(edgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
    EI(edgeIdx, 0) = IGL_COLLAPSE_EDGE_NULL;
    EI(edgeIdx, 1) = IGL_COLLAPSE_EDGE_NULL;
  };

  // update edge info for each flap
  const int m = tris.rows();
  for(int side = 0;side<2;side++)
  {
    const int f = EF(edgeIdx, side);
    const int v = EI(edgeIdx, side);
    const int sign = (eflip==0?1:-1)*(1-2*side);
    // next edge emanating from d
    const int e1 = EMAP(f+m*((v+sign*1+3)%3));
    // prev edge pointing to s
    const int e2 = EMAP(f+m*((v+sign*2+3)%3));
    assert(edges(e1, 0) == d || edges(e1, 1) == d);
    assert(edges(e2, 0) == s || edges(e2, 1) == s);
    // face adjacent to f on e1,  also incident on d
    const bool flip1 = EF(e1, 1)==f;
    const int f1 = flip1 ? EF(e1, 0) : EF(e1, 1);
    assert(f1!=f);
    assert(tris(f1, 0)==d || tris(f1, 1)==d || tris(f1, 2) == d);
    // across from which vertex of f1 does e1 appear?
    const int v1 = flip1 ? EI(e1, 0) : EI(e1, 1);
    // Kill e1
    kill_edge(e1);
    // Kill f
    tris(f, 0) = IGL_COLLAPSE_EDGE_NULL;
    tris(f, 1) = IGL_COLLAPSE_EDGE_NULL;
    tris(f, 2) = IGL_COLLAPSE_EDGE_NULL;
    // map f1's edge on e1 to e2
    assert(EMAP(f1+m*v1) == e1);
    EMAP(f1+m*v1) = e2;
    // side opposite f2,  the face adjacent to f on e2,  also incident on s
    const int opp2 = (EF(e2, 0)==f?0:1);
    assert(EF(e2, opp2) == f);
    EF(e2, opp2) = f1;
    EI(e2, opp2) = v1;
    // remap e2 from d to s
    edges(e2, 0) = edges(e2, 0)==d ? s : edges(e2, 0);
    edges(e2, 1) = edges(e2, 1)==d ? s : edges(e2, 1);
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

  // finally,  reindex faces and edges incident on d. Do this last so asserts make sense.
  // Could actually skip first and last,  since those are always the two
  // collpased faces. Nah,  this is handled by (tris(f, v) == d)
  // Don't attempt to use Nde, Nse here because EMAP has changed
  {
    int p1 = -1;
    for(auto f : nV2Fd)
    {
      for(int v = 0;v<3;v++)
      {
        if(tris(f, v) == d)
        {
          const int e1 = EMAP(f+m*((v+1)%3));
          const int flip1 = (EF(e1, 0)==f)?1:0;
          assert( edges(e1, flip1) == d || edges(e1, flip1) == s);
          edges(e1, flip1) = s;
          const int e2 = EMAP(f+m*((v+2)%3));
          // Skip if we just handled this edge (claim: this will be all except for the first non-trivial face)
          if(e2 != p1)
          {
            const int flip2 = (EF(e2, 0)==f)?0:1;
            assert( edges(e2, flip2) == d || edges(e2, flip2) == s);
            edges(e2, flip2) = s;
          }

          tris(f, v) = s;
          p1 = e1;
          break;
        }
      }
    }
  }
  // Finally,  "remove" this edge and its information
  kill_edge(edgeIdx);
  return true;
}


// 重载2.1
IGL_INLINE bool igl::collapse_edge(
  const int edgeIdx, 
  const Eigen::RowVectorXd & p, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI)
{
  int e1, e2, f1, f2;
  return collapse_edge(edgeIdx, p, vers, tris, edges, EMAP, EF, EI, e1, e2, f1, f2);
}


// 重载2.2
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & EQ, 
  Eigen::MatrixXd & C)
{
  int edgeIdx, e1, e2, f1, f2;
  decimate_pre_collapse_callback always_try;
  decimate_post_collapse_callback never_care;
  decimate_trivial_callbacks(always_try, never_care);
  return 
    collapse_edge(cost_and_placement, always_try, never_care, \
            vers, tris, edges, EMAP, EF, EI, workQueue, EQ, C, edgeIdx, e1, e2, f1, f2);
}


// 重载2.3
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_pre_collapse_callback       & pre_collapse, 
  const decimate_post_collapse_callback      & post_collapse, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & EQ, 
  Eigen::MatrixXd & C)
{
  int edgeIdx, e1, e2, f1, f2;
  return 
    collapse_edge(
      cost_and_placement, pre_collapse, post_collapse, 
      vers, tris, edges, EMAP, EF, EI, workQueue, EQ, C, edgeIdx, e1, e2, f1, f2);
}


// 重载2.
IGL_INLINE bool igl::collapse_edge(
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_pre_collapse_callback       & pre_collapse, 
  const decimate_post_collapse_callback      & post_collapse, 
  Eigen::MatrixXd & vers, 
  Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & edges, 
  Eigen::VectorXi & EMAP, 
  Eigen::MatrixXi & EF, 
  Eigen::MatrixXi & EI, 
  igl::min_heap< std::tuple<double, int, int> > & workQueue, 
  Eigen::VectorXi & EQ, 
  Eigen::MatrixXd & C, 
  int & edgeIdx, 
  int & e1, 
  int & e2, 
  int & f1, 
  int & f2)
{
  using namespace Eigen;
  using namespace igl;

  std::tuple<double, int, int> p;
  while(true)
  {
    if(workQueue.empty())
    {
      // no edges to collapse
      edgeIdx = -1;
      return false;
    }

    // pop from workQueue
    p = workQueue.top();
    if(std::get<0>(p) == std::numeric_limits<double>::infinity())
    {
      edgeIdx = -1;     // min cost edge is infinite cost
      return false;
    }
    workQueue.pop();
    edgeIdx = std::get<1>(p);

    // Check if matches timestamp
    if(std::get<2>(p) == EQ(edgeIdx))
         break;

    assert(std::get<2>(p)  < EQ(edgeIdx) || EQ(edgeIdx) == -1);         // must be stale or dead.

    // try again.
  }


  // Why is this computed up here?
  
  // If we just need original face neighbors of edge,  could we gather that more directly than gathering face neighbors of each vertex?
  std::vector<int> /*Nse, */Nsf, Nsv;
  circulation(edgeIdx,  true, tris, EMAP, EF, EI, /*Nse, */Nsv, Nsf);
  std::vector<int> /*Nde, */Ndf, Ndv;
  circulation(edgeIdx,  false, tris, EMAP, EF, EI, /*Nde, */Ndv, Ndf);


  bool collapsed = true;
  if(pre_collapse(vers, tris, edges, EMAP, EF, EI, workQueue, EQ, C, edgeIdx))
  {
    collapsed = collapse_edge(
      edgeIdx, C.row(edgeIdx), 
      Nsv, Nsf, Ndv, Ndf, 
      vers, tris, edges, EMAP, EF, EI, e1, e2, f1, f2);
  }
  else
     collapsed = false;          // Aborted by pre collapse callback


  post_collapse(vers, tris, edges, EMAP, EF, EI, workQueue, EQ, C, edgeIdx, e1, e2, f1, f2, collapsed);
  if(collapsed)
  {
    // Erase the two,  other collapsed edges by marking their timestamps as -1
    EQ(e1) = -1;
    EQ(e2) = -1;
    // TODO: visits edges multiple times,  ~150% more updates than should
    //
    // update local neighbors
    // loop over original face neighbors
    //
    // Can't use previous computed Nse and Nde because those refer to EMAP
    // before it was changed...
    std::vector<int> Nf;
    Nf.reserve( Nsf.size() + Ndf.size() ); // preallocate memory
    Nf.insert( Nf.end(),  Nsf.begin(),  Nsf.end() );
    Nf.insert( Nf.end(),  Ndf.begin(),  Ndf.end() );
    // https://stackoverflow.com/a/1041939/148668
    std::sort( Nf.begin(),  Nf.end() );
    Nf.erase( std::unique( Nf.begin(),  Nf.end() ),  Nf.end() );
    // Collect all edges that must be updated
    std::vector<int> Ne;
    Ne.reserve(3*Nf.size());
    for(auto & n : Nf)
    {
      if(tris(n, 0) != IGL_COLLAPSE_EDGE_NULL ||
          tris(n, 1) != IGL_COLLAPSE_EDGE_NULL ||
          tris(n, 2) != IGL_COLLAPSE_EDGE_NULL)
      {
        for(int v = 0;v<3;v++)
        {
          // get edge id
          const int ei = EMAP(v*tris.rows()+n);
          Ne.push_back(ei);
        }
      }
    }
    // Only process edge once
    std::sort( Ne.begin(),  Ne.end() );
    Ne.erase( std::unique( Ne.begin(),  Ne.end() ),  Ne.end() );
    for(auto & ei : Ne)
    {
       // compute cost and potential placement
       double cost;
       RowVectorXd place;
       cost_and_placement(ei, vers, tris, edges, EMAP, EF, EI, cost, place);
       // Increment timestamp
       EQ(ei)++;
       // Replace in queue
       workQueue.emplace(cost, ei, EQ(ei));
       C.row(ei) = place;
    }
  }
  else
  {
    // reinsert with infinite weight (the provided cost function must **not**
    // have given this un-collapsable edge inf cost already)
    // Increment timestamp
    EQ(edgeIdx)++;
    // Replace in queue
    workQueue.emplace(std::numeric_limits<double>::infinity(), edgeIdx, EQ(edgeIdx));
  }


  return collapsed;
}
