#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/exact_geodesic.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/parula.h>
#include <igl/isolines_map.h>
#include <igl/PI.h>
#include <iostream>
#include "tutorial_shared_path.h"


// 测量测地距离；

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  Eigen::MatrixXd vers;
  Eigen::MatrixXi F;
  igl::opengl::glfw::Viewer viewer;

  igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", vers, F);

  const auto update_distance = [&](const int vid)
  {
    Eigen::VectorXi VS,FS,VT,FT;

    // The selected vertex is the source
    VS.resize(1);
    VS << vid;

    // All vertices are the targets
    VT.setLinSpaced(vers.rows(),0,vers.rows()-1);
    Eigen::VectorXd d;
    std::cout<<"Computing geodesic distance to vertex "<<vid<<"..."<<std::endl;
    igl::exact_geodesic(vers,F,VS,FS,VT,FT,d);

    // Plot the mesh
    Eigen::MatrixXd CM;
    igl::parula(Eigen::VectorXd::LinSpaced(21,0,1).eval(),false,CM);
    igl::isolines_map(Eigen::MatrixXd(CM),CM);
    viewer.data().set_colormap(CM);
    viewer.data().set_data(d);

    std::cout << "Computing finished." << std::endl;
  };

  // 鼠标左键点击一个顶点，计算网格上其他顶点到该顶点的测地距离，用不同色彩表示距离的远近。
  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(
      Eigen::Vector2f(x,y),
      viewer.core().view,
      viewer.core().proj,
      viewer.core().viewport,
      vers,
      F,
      fid,
      bc))
    {
      int max;
      bc.maxCoeff(&max);
      int vid = F(fid,max);
      update_distance(vid);
      return true;
    }
    return false;
  };

  viewer.data().set_mesh(vers,F);
  viewer.data().show_lines = false;

  cout << "Click on mesh to define new source.\n" << std::endl;

  update_distance(0);
  return viewer.launch();
}
