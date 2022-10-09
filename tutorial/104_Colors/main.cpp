#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd vers;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/screwdriver.off", vers, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, F);

  // Use the (normalized) vertex positions as colors
  C = 
    (vers.rowwise()            - vers.colwise().minCoeff()).array().rowwise()/
    (vers.colwise().maxCoeff() - vers.colwise().minCoeff()).array();

  // Add per-vertex colors
  viewer.data().set_colors(C);

  // Launch the viewer
  viewer.launch();
}
