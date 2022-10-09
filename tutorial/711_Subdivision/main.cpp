#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi OF,F;
  Eigen::MatrixXd OV,vers;
  bool show_swept_volume = false;
  read_triangle_mesh(
      TUTORIAL_SHARED_PATH "/decimated-knight.off",OV,OF);
  vers = OV;
  F = OF;
  cout<<R"(Usage:
1  Restore Original mesh
2  Apply In-plane upsampled mesh
3  Apply Loop subdivided mesh
4  Apply False barycentric subdivision
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers,F);
  viewer.data().set_face_based(true);

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
        {
          vers = OV;
          F = OF;
          break;
        }
        case '2':
        {
          igl::upsample( Eigen::MatrixXd(vers), Eigen::MatrixXi(F), vers,F);
          break;
        }
        case '3':
        {
          igl::loop( Eigen::MatrixXd(vers), Eigen::MatrixXi(F), vers,F);
          break;
        }
        case '4':
        {
          igl::false_barycentric_subdivision(
            Eigen::MatrixXd(vers),Eigen::MatrixXi(F),vers,F);
          break;
        }
      }
      viewer.data().clear();
      viewer.data().set_mesh(vers,F);
      viewer.data().set_face_based(true);
      return true;
    };
  viewer.launch();
}
