#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd vers;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", vers, F);

  igl::opengl::glfw::Viewer viewer;			// libigl�л���glfw��ͼ�ο�ܡ�
  viewer.data().set_mesh(vers, F);

  // Ĭ����������ת������������һ���趨Ϊ������ת��
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}
