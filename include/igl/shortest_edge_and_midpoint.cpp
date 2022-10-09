#include "shortest_edge_and_midpoint.h"

IGL_INLINE void igl::shortest_edge_and_midpoint(
  const int edgeIdx,
  const Eigen::MatrixXd & vers,
  const Eigen::MatrixXi & /*tris*/,
  const Eigen::MatrixXi & uEdges,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & edgeCenter)
{
	cost = (vers.row(uEdges(edgeIdx,0))-vers.row(uEdges(edgeIdx,1))).norm();
	edgeCenter = 0.5*(vers.row(uEdges(edgeIdx,0))+vers.row(uEdges(edgeIdx,1)));
}
