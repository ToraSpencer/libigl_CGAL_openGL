#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  using namespace Eigen;
  std::string filename = TUTORIAL_SHARED_PATH "/fertility.off";
  if(argc>1)
  {
    filename = argv[1];
  }


  igl::read_triangle_mesh(filename, V, F);

  // 1. 分布计算平均曲率：
  MatrixXd HN;
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);
  // Laplace-Beltrami of position
  HN = -Minv*(L*V);
  // Extract magnitude as mean curvature
  VectorXd H1 = HN.rowwise().norm();


  // 2. 直接计算——igl::principal_curvature()——使用二次曲面拟合计算最大、最小曲率及其方向
  MatrixXd PD1,PD2;
  VectorXd PV1,PV2;
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);

  VectorXd H2 = 0.5*(PV1+PV2);          // 平均曲率

  // 3. 比较分布计算和直接计算的平均曲率结果：
  std::cout << "(H1 - H2).norm() == " << (H1 - H2).norm() << std::endl;

  // 4. 绘图：
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(H2);
  const double avg = igl::avg_edge_length(V,F);         // 所有边的平均边长，作为指示线的长度
  const RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
  viewer.data().add_edges(V + PD1*avg, V - PD1*avg, red);         // 最大曲率方向用红色指示线标识
  viewer.data().add_edges(V + PD2*avg, V - PD2*avg, blue);        // 最小曲率方向用蓝色指示线标识

  viewer.data().show_lines = false;                 // 隐藏网格线

  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // 设定可以三轴旋转
  viewer.launch();
}
