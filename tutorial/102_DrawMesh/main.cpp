#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd vers;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", vers, F);

  igl::opengl::glfw::Viewer viewer;			// libigl中基于glfw的图形框架。
  viewer.data().set_mesh(vers, F);

  // 默认是两轴旋转，加上下面这一句设定为三轴旋转。
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}
