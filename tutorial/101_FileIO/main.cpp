#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;


int main(int argc, char *argv[])
{
	std::cout << "hello world" << std::endl;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH  "/cube.off", V, F);


  // Print the vertices and faces matrices
  std::cout << "Vertices: " << std::endl << V << std::endl;
  std::cout << "Faces:    " << std::endl << F << std::endl;


  // Save the mesh in OBJ format
  igl::writeOBJ("cube.obj",V, F);

  igl::readOBJ("cube.obj", V, F);
}
