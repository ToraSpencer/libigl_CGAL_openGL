#include "qslim.h"

#include "collapse_edge.h"
#include "connect_boundary_to_infinity.h"
#include "decimate.h"
#include "edge_flaps.h"
#include "is_edge_manifold.h"
#include "max_faces_stopping_condition.h"
#include "per_vertex_point_to_plane_quadrics.h"
#include "qslim_optimal_collapse_edge_callbacks.h"
#include "quadric_binary_plus_operator.h"
#include "remove_unreferenced.h"
#include "slice.h"
#include "slice_mask.h"


IGL_INLINE bool igl::qslim(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const size_t max_m, 
  Eigen::MatrixXd & versOut, 
  Eigen::MatrixXi & trisOut, 
  Eigen::VectorXi & newOldTrisInfo, 
  Eigen::VectorXi & newOldVersInfo)
{
  using namespace igl; 
 
  const int trisCountOri = tris.rows(); 
  int trisCount = tris.rows(); 
  typedef Eigen::MatrixXd DerivedV; 
  typedef Eigen::MatrixXi DerivedF; 
  DerivedV vers0; 
  DerivedF tris0; 
  Eigen::VectorXi edgeUeInfo;
  Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
  decimate_cost_and_placement_callback cost_and_placement;
  decimate_pre_collapse_callback       pre_collapse;
  decimate_post_collapse_callback      post_collapse;
  typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
  std::vector<Quadric> quadrics;
  int v1 = -1;       // State variables keeping track of edge we just collapsed
  int v2 = -1;

  // 1. 
  igl::connect_boundary_to_infinity(vers, tris,  vers0,  tris0); 

  // this includes meshes with non-manifold vertices on the boundary since these will create a non-manifold edge when connected to infinity.
  if(!is_edge_manifold(tris0))       // decimate will not work correctly on non-edge-manifold meshes. 
        return false; 

  edge_flaps(tris0, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo); 


  // 2. 计算每个顶点的Q矩阵：
  per_vertex_point_to_plane_quadrics(vers0, tris0, edgeUeInfo, UeTrisInfo, UeCornersInfo, quadrics); 


  // 3. 获取边折叠前后的pre_collapse(), post_collapse()函数子，以及计算边折叠cost值的函数子cost_and_placement();
  qslim_optimal_collapse_edge_callbacks(uEdges, quadrics, v1, v2,  cost_and_placement, pre_collapse, post_collapse); 
  
  // 4. Call to greedy decimator
  bool ret = decimate(vers0,  tris0, 
        cost_and_placement, 
        max_faces_stopping_condition(trisCount, trisCountOri, max_m), 
        pre_collapse, 
        post_collapse, 
        uEdges,  edgeUeInfo,  UeTrisInfo,  UeCornersInfo, 
        versOut,  trisOut,  newOldTrisInfo,  newOldVersInfo); 

  // 5. Remove phony boundary faces and clean up
  const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCountOri); 
  igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut); 
  igl::slice_mask(Eigen::VectorXi(newOldTrisInfo),  keep,  1,  newOldTrisInfo); 
  Eigen::VectorXi _1, I2; 
  igl::remove_unreferenced(Eigen::MatrixXd(versOut),  Eigen::MatrixXi(trisOut),  versOut,  trisOut, _1, I2); 
  igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo); 

  return ret; 
}
