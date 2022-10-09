#include "edge_flaps.h"
#include "unique_edge_map.h"
#include <vector>
#include <cassert>


// жиди1
IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & tris, 
  const Eigen::MatrixXi & uEdges, 
  const Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & UeTrisInfo, 
  Eigen::MatrixXi & UeCornersInfo)
{
  // Initialize to boundary value
  UeTrisInfo.setConstant(uEdges.rows(), 2, -1);
  UeCornersInfo.setConstant(uEdges.rows(), 2, -1);

  // loop over all faces
  for(int f = 0;f<tris.rows();f++)
  {
    // loop over edges across from corners
    for(int v = 0;v<3;v++)
    {
      // get edge id
      const int e = edgeUeInfo(v*tris.rows()+f);
      // See if this is left or right flap w.r.t. edge orientation
      if( tris(f, (v+1)%3) == uEdges(e, 0) && tris(f, (v+2)%3) == uEdges(e, 1))
      {
        UeTrisInfo(e, 0) = f;
        UeCornersInfo(e, 0) = v;
      }else
      {
        assert(tris(f, (v+1)%3) == uEdges(e, 1) && tris(f, (v+2)%3) == uEdges(e, 0));
        UeTrisInfo(e, 1) = f;
        UeCornersInfo(e, 1) = v;
      }
    }
  }
}


// жиди1.1
IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & tris, 
  Eigen::MatrixXi & uEdges, 
  Eigen::VectorXi & edgeUeInfo, 
  Eigen::MatrixXi & UeTrisInfo, 
  Eigen::MatrixXi & UeCornersInfo)
{
  Eigen::MatrixXi allE;
  igl::unique_edge_map(tris, allE, uEdges, edgeUeInfo);
  // Const-ify to call overload
  const auto & cuE = uEdges;
  const auto & cEMAP = edgeUeInfo;
  return edge_flaps(tris, cuE, cEMAP, UeTrisInfo, UeCornersInfo);
}
