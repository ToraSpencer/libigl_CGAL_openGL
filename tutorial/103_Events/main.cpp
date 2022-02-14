#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V1,V2;
Eigen::MatrixXi F1,F2;


// 键盘事件回调函数——切换显示两个不同的网格：
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
      if (key == '1')
      {
        viewer.data().clear();
        viewer.data().set_mesh(V1, F1);
        viewer.core().align_camera_center(V1,F1);
      }
      else if (key == '2')
      {
        viewer.data().clear();
        viewer.data().set_mesh(V2, F2);
        viewer.core().align_camera_center(V2,F2);
      }

  return false;
}
 

int main(int argc, char *argv[])
{
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V1, F1);
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V2, F2);
  std::cout << R"(
    1 Switch to bump mesh
    2 Switch to fertility mesh
    )";

  igl::opengl::glfw::Viewer viewer;

  viewer.callback_key_down = &key_down;          
  viewer.data().set_mesh(V1, F1);
  viewer.launch();
}
