#include <igl/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MatrixXi tris;
  MatrixXd vers;

  read_triangle_mesh(TUTORIAL_SHARED_PATH "/armadillo.obj",vers,tris);
  cout<<"Creating grid..."<<endl;

  // number of vertices on the largest side
  const int s = 100;

  // 生成栅格：
  MatrixXd GV;
  Eigen::RowVector3i res;
  igl::voxel_grid(vers,0,s,1,GV,res);
 
  // compute values
  cout<<"Computing distances..."<<endl;

  VectorXd S,B;

  {
    VectorXi I;
    MatrixXd C,N;
    signed_distance(GV,vers,tris,SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);

    // 符号距离场改写为符号场——网格内为-1，网格面上为0，外面为1：
    B = S;
    for_each(B.data(),B.data()+B.size(),[](double& b)\
        {
            b=(b>0 ? 1 : (b<0 ? -1 : 0));
        });
  }

  MatrixXd SV,BV;
  MatrixXi SF,BF;
  igl::marching_cubes(S,GV,res(0),res(1),res(2),0,SV,SF);
  igl::marching_cubes(B,GV,res(0),res(1),res(2),0,BV,BF);

  cout<<R"(Usage:
'1'  Show original mesh.
'2'  Show marching cubes contour of signed distance.
'3'  Show marching cubes contour of indicator function.
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(SV,SF);
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
          viewer.data().clear();
          viewer.data().set_mesh(vers,tris);
          break;
        case '2':
          viewer.data().clear();
          viewer.data().set_mesh(SV,SF);
          break;
        case '3':
          viewer.data().clear();
          viewer.data().set_mesh(BV,BF);
          break;
      }
      viewer.data().set_face_based(true);
      return true;
    };
  viewer.launch();
}
