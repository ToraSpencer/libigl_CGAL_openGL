#include "edge_collapse_is_valid.h"
#include "collapse_edge.h"
#include "circulation.h"
#include "intersect.h"
#include "unique.h"
#include "list_to_matrix.h"
#include <vector>


// жиди1
IGL_INLINE bool igl::edge_collapse_is_valid(
  const int e,
  const Eigen::MatrixXi & tris,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & edgeUeInfo,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI)
{
  using namespace Eigen;
  using namespace std;
  // For consistency with collapse_edge.cpp, let's determine edge flipness
  // (though not needed to check validity)
  const int eflip = E(e,0)>E(e,1);
  // source and destination
  const int s = eflip?E(e,1):E(e,0);
  const int d = eflip?E(e,0):E(e,1);

  if(s == IGL_COLLAPSE_EDGE_NULL && d==IGL_COLLAPSE_EDGE_NULL)
  {
    return false;
  }
  // check if edge collapse is valid: intersection of vertex neighbors of s and
  // d should be exactly 2+(s,d) = 4
  // http://stackoverflow.com/a/27049418/148668
  {
    // all vertex neighbors around edge, including the two vertices of the edge
    const auto neighbors = [](
      const int e,
      const bool ccw,
      const Eigen::MatrixXi & tris,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & edgeUeInfo,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI) 
    {
      vector<int> N,uN;
      vector<int> V2Fe = circulation(e, ccw,edgeUeInfo,EF,EI);
      for(auto f : V2Fe)
      {
        N.push_back(tris(f,0));
        N.push_back(tris(f,1));
        N.push_back(tris(f,2));
      }
      vector<size_t> _1,_2;
      igl::unique(N,uN,_1,_2);
      VectorXi uNm;
      list_to_matrix(uN,uNm);
      return uNm;
    };
    VectorXi Ns = neighbors(e, eflip,tris,E,edgeUeInfo,EF,EI);
    VectorXi Nd = neighbors(e,!eflip,tris,E,edgeUeInfo,EF,EI);
    VectorXi Nint = igl::intersect(Ns,Nd);
    if(Nint.size() != 4)
    {
      return false;
    }
    if(Ns.size() == 4 && Nd.size() == 4)
    {
      VectorXi NsNd(8);
      NsNd<<Ns,Nd;
      VectorXi Nun,_1,_2;
      igl::unique(NsNd,Nun,_1,_2);
      // single tet, don't collapse
      if(Nun.size() == 4)
      {
        return false;
      }
    }
  }
  return true;
}


// жиди2ЃЛ
IGL_INLINE bool igl::edge_collapse_is_valid(
  std::vector<int> & srcNbrIdx,
  std::vector<int> & desNbrIdx)
{
  // Do we really need to check if edge is IGL_COLLAPSE_EDGE_NULL ?
  if(srcNbrIdx.size()<2 || desNbrIdx.size()<2)
  {
    // Bogus data
    assert(false);
    return false;
  }

  // determine if the first two vertices are the same before reordering.
  // If they are and there are 3 each, then (I claim) this is an edge on a single tet.
  const bool first_two_same = (srcNbrIdx[0] == desNbrIdx[0]) && (srcNbrIdx[1] == desNbrIdx[1]);
  if(srcNbrIdx.size() == 3 && desNbrIdx.size() == 3 && first_two_same) 
    return false;           // single tet
 

  // https://stackoverflow.com/a/19483741/148668
  std::sort(srcNbrIdx.begin(), srcNbrIdx.end());
  std::sort(desNbrIdx.begin(), desNbrIdx.end());
  std::vector<int> Nint;
  std::set_intersection(srcNbrIdx.begin(), srcNbrIdx.end(), desNbrIdx.begin(), desNbrIdx.end(), std::back_inserter(Nint));

  // check if edge collapse is valid: intersection of vertex neighbors of s and d should be exactly 2+(s,d) = 4

  // http://stackoverflow.com/a/27049418/148668
  if(Nint.size() != 2)
    return false;
 
  
  return true;
}
