#include "edge_flaps.h"
#include "unique_edge_map.h"
#include <vector>
#include <cassert>

IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & tris,
  const Eigen::MatrixXi & uE,
  const Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  // Initialize to boundary value
  EF.setConstant(uE.rows(),2,-1);
  EI.setConstant(uE.rows(),2,-1);
  // loop over all faces
  for(int f = 0;f<tris.rows();f++)
  {
    // loop over edges across from corners
    for(int v = 0;v<3;v++)
    {
      // get edge id
      const int e = EMAP(v*tris.rows()+f);
      // See if this is left or right flap w.r.t. edge orientation
      if( tris(f,(v+1)%3) == uE(e,0) && tris(f,(v+2)%3) == uE(e,1))
      {
        EF(e,0) = f;
        EI(e,0) = v;
      }else
      {
        assert(tris(f,(v+1)%3) == uE(e,1) && tris(f,(v+2)%3) == uE(e,0));
        EF(e,1) = f;
        EI(e,1) = v;
      }
    }
  }
}

IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & tris,
  Eigen::MatrixXi & uE,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  Eigen::MatrixXi allE;
  igl::unique_edge_map(tris,allE,uE,EMAP);
  // Const-ify to call overload
  const auto & cuE = uE;
  const auto & cEMAP = EMAP;
  return edge_flaps(tris,cuE,cEMAP,EF,EI);
}
