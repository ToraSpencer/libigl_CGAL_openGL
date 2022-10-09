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

Eigen::MatrixXd vers;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  using namespace Eigen;

  std::string filename = TUTORIAL_SHARED_PATH "/fertility.off";
  if(argc>1)
    filename = argv[1];

  igl::read_triangle_mesh(filename, vers, F);

  // 1. 手动分步计算平均曲率：
  MatrixXd HN;
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(vers,F,L);
  igl::massmatrix(vers,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);

  //    Laplace-Beltrami of position
  HN = -Minv*(L*vers);

  //    Extract magnitude as mean curvature
  VectorXd H1 = HN.rowwise().norm();

  // 2. 直接使用轮子计算 
  MatrixXd PD1,PD2;
  VectorXd PV1,PV2;

  //    igl::principal_curvature()——使用二次曲面拟合计算最大、最小曲率及其方向
  igl::principal_curvature(vers,F,PD1,PD2,PV1,PV2);

  VectorXd H2 = 0.5*(PV1+PV2);          // 平均曲率

  // 3. 比较手动分步计算和直接计算的平均曲率结果：
  std::cout << "(H1 - H2).norm() == " << (H1 - H2).norm() << std::endl;

  // 4. 绘图：
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, F);
  viewer.data().set_data(H2);

  const double aveLen = igl::avg_edge_length(vers,F);         // 所有边的平均边长，作为指示线的长度
  const RowVector3d red(0.8,0.2,0.2), blue(0.2,0.2,0.8);    // RGB色彩向量；

  //        指示线用边数据的形式渲染出来；
  viewer.data().add_edges(vers + PD1*aveLen, vers - PD1*aveLen, red);         // 最大曲率方向用红色指示线标识
  viewer.data().add_edges(vers + PD2*aveLen, vers - PD2*aveLen, blue);        // 最小曲率方向用蓝色指示线标识

  viewer.data().show_lines = false;                 // 隐藏网格线

  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // 设定可以三轴旋转
  viewer.launch();
}
