#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd vers;
Eigen::MatrixXi F;

Eigen::MatrixXd N_vertices;
Eigen::MatrixXd N_faces;
Eigen::MatrixXd N_corners;


// 键盘事件回调函数
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  switch(key)
  {
    case '1':
      viewer.data().set_normals(N_faces);
      return true;
    case '2':
      viewer.data().set_normals(N_vertices);
      return true;
    case '3':
      viewer.data().set_normals(N_corners);
      return true;
    default: break;
  }
  return false;
}


int main(int argc, char *argv[])
{
  igl::read_triangle_mesh(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/fandisk.off",vers,F);

  // 面法线向量
  igl::per_face_normals(vers,F,N_faces);

  // 顶点法向
  igl::per_vertex_normals(vers,F,N_vertices);

  // Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
  igl::per_corner_normals(vers, F, 20, N_corners);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = true;
  viewer.data().set_mesh(vers, F);
  viewer.data().set_normals(N_faces);

  std::cout<<
    "Press '1' for per-face normals."<<std::endl<<
    "Press '2' for per-vertex normals."<<std::endl<<
    "Press '3' for per-corner normals."<<std::endl;
  viewer.launch();
}
