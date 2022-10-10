#include "circulation.h"
#include "list_to_matrix.h"
#include <cassert>

// 重载1
IGL_INLINE std::vector<int> igl::circulation(
  const int edgeIdx, 
  const bool ccw, 
  const Eigen::VectorXi & EMAP, 
  const Eigen::MatrixXi & EF, 
  const Eigen::MatrixXi & EI)
{
  // prepare output
  std::vector<int> N;
  N.reserve(6);
  const int m = EMAP.size()/3;
  assert(m*3 == EMAP.size());
  const auto & step = [&](const int edgeIdx,  const int ff, int & ne, int & nf)
  {
    assert((EF(edgeIdx, 1) == ff || EF(edgeIdx, 0) == ff) && "edgeIdx should touch ff");
    //const int fside = EF(edgeIdx, 1)==ff?1:0;
    const int nside = EF(edgeIdx, 0)==ff?1:0;
    const int nv = EI(edgeIdx, nside);
    // get next face
    nf = EF(edgeIdx, nside);
    // get next edge 
    const int dir = ccw?-1:1;
    ne = EMAP(nf+m*((nv+dir+3)%3));
  };
  // Always start with first face (ccw in step will be sure to turn right
  // direction)
  const int f0 = EF(edgeIdx, 0);
  int fi = f0;
  int ei = edgeIdx;
  while(true)
  {
    step(ei, fi, ei, fi);
    N.push_back(fi);
    // back to start?
    if(fi == f0)
    {
      assert(ei == edgeIdx);
      break;
    }
  }
  return N;
}


// 重载1.1
IGL_INLINE void igl::circulation(
  const int edgeIdx, 
  const bool ccw, 
  const Eigen::VectorXi & EMAP, 
  const Eigen::MatrixXi & EF, 
  const Eigen::MatrixXi & EI, 
  Eigen::VectorXi & vN)
{
  std::vector<int> N = circulation(edgeIdx, ccw, EMAP, EF, EI);
  igl::list_to_matrix(N, vN);
}


// 重载2；
IGL_INLINE void igl::circulation(
  const int edgeIdx, 
  const bool ccw, 
  const Eigen::MatrixXi & tris, 
  const Eigen::VectorXi & EMAP, 
  const Eigen::MatrixXi & EF, 
  const Eigen::MatrixXi & EI, 
  std::vector<int> & nbrVersIdx, 
  std::vector<int> & nbrTrisIdx)
{
  /*
       for edgeIdx --> (bf) and ccw=true
  
           c---d
          / \ / \
         a---b-edgeIdx-f
          \ / \ /
           g---h
  
        // (might start with {bhf} depending on edge)
        Ne = […] -> [fd db dc cb ca ab ag gb gh hb hf fb]
                    {upto cylic order}
        nbrTrisIdx = […] -> [{bfd},  {bdc},  {bca},  {bag},  {bgh},  {bhf}]
        nbrVersIdx = [d c a g h f]
  
       prepare output
      Ne.clear();Ne.reserve(2*10);
    */

  nbrVersIdx.clear(); nbrVersIdx.reserve(10);
  nbrTrisIdx.clear(); nbrTrisIdx.reserve(10);
  const int m = EMAP.size()/3;
  assert(m*3 == EMAP.size());

  const auto & step = [&]( const int edgeIdx,  const int ff,   int & ne,   /*int& re,  */  int & rv,   int & nf)
  {
    assert((EF(edgeIdx, 1) == ff || EF(edgeIdx, 0) == ff) && "edgeIdx should touch ff");
    const int nside = EF(edgeIdx, 0)==ff?1:0;
    const int nv = EI(edgeIdx, nside);

    // get next face
    nf = EF(edgeIdx, nside);

    // get next edge 
    const int dir = ccw?-1:1;
    rv = tris(nf,  nv);
    ne = EMAP(nf+m*((nv+dir+3)%3));
    //re = EMAP(nf+m*((nv+2*dir+3)%3));
  };

  // Always start with first face (ccw in step will be sure to turn right direction)
  const int f0 = EF(edgeIdx, 0);
  int fi = f0;
  int ei = edgeIdx;

  while(true)
  {
    int re, rv;
    step(ei, fi, ei/*, re*/, rv, fi);
    nbrTrisIdx.push_back(fi);
    nbrVersIdx.push_back(rv);

    // back to start?
    if(fi == f0)
    {
      assert(ei == edgeIdx);
      break;
    }
  }
}
