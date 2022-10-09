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

  // 1. �ֶ��ֲ�����ƽ�����ʣ�
  MatrixXd HN;
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(vers,F,L);
  igl::massmatrix(vers,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);

  //    Laplace-Beltrami of position
  HN = -Minv*(L*vers);

  //    Extract magnitude as mean curvature
  VectorXd H1 = HN.rowwise().norm();

  // 2. ֱ��ʹ�����Ӽ��� 
  MatrixXd PD1,PD2;
  VectorXd PV1,PV2;

  //    igl::principal_curvature()����ʹ�ö���������ϼ��������С���ʼ��䷽��
  igl::principal_curvature(vers,F,PD1,PD2,PV1,PV2);

  VectorXd H2 = 0.5*(PV1+PV2);          // ƽ������

  // 3. �Ƚ��ֶ��ֲ������ֱ�Ӽ����ƽ�����ʽ����
  std::cout << "(H1 - H2).norm() == " << (H1 - H2).norm() << std::endl;

  // 4. ��ͼ��
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, F);
  viewer.data().set_data(H2);

  const double aveLen = igl::avg_edge_length(vers,F);         // ���бߵ�ƽ���߳�����Ϊָʾ�ߵĳ���
  const RowVector3d red(0.8,0.2,0.2), blue(0.2,0.2,0.8);    // RGBɫ��������

  //        ָʾ���ñ����ݵ���ʽ��Ⱦ������
  viewer.data().add_edges(vers + PD1*aveLen, vers - PD1*aveLen, red);         // ������ʷ����ú�ɫָʾ�߱�ʶ
  viewer.data().add_edges(vers + PD2*aveLen, vers - PD2*aveLen, blue);        // ��С���ʷ�������ɫָʾ�߱�ʶ

  viewer.data().show_lines = false;                 // ����������

  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);    // �趨����������ת
  viewer.launch();
}
