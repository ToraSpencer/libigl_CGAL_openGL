#include <igl/barycenter.h>
#include <igl/colon.h>
#include <igl/jet.h>
#include <igl/readOFF.h>
#include <igl/slice_into.h>
#include <igl/sortrows.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd vers;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off",vers,F);

  // Sort barycenters lexicographically
  MatrixXd BC,sorted_BC;
  igl::barycenter(vers,F,BC);
  VectorXi I,J;
  // sorted_BC = BC(I,:)
  igl::sortrows(BC,true,sorted_BC,I);
  // Get sorted "place" from sorted indices
  J.resize(I.rows());
  // J(I) = 1:numel(I)
  igl::slice_into(igl::colon<int>(0,I.size()-1),I,J);

  // Pseudo-color based on sorted place
  MatrixXd C;
  igl::jet(J,true,C);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, F);
  viewer.data().set_colors(C);
  viewer.launch();
}
