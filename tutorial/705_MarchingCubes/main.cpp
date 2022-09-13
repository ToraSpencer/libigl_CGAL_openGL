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

  read_triangle_mesh(TUTORIAL_SHARED_PATH "/armadillo.obj", vers, tris);
  cout<<"Creating grid..."<<endl;

  // const int num = 100;              // number of vertices on the largest side
  const int num = 10;              // number of vertices on the largest side

  // 1. ����դ��
  MatrixXd gridCenters;
  Eigen::RowVector3i gridCounts;
  igl::voxel_grid(vers, 0, num, 1, gridCenters, gridCounts);
  
  cout<<"Computing distances..."<<endl;           // compute values

  VectorXd SDF, signValues;

  {
    VectorXi I;                     // useless
    MatrixXd C, N;              // useless

    // 2. ������ž��볡
    signed_distance(gridCenters, vers, tris, SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N); 

    // 3. ���ž��볡��дΪ���ų�����������Ϊ-1����������Ϊ0������Ϊ1��
    signValues = SDF;
    for_each(signValues.data(), signValues.data()+signValues.size(), [](double& b)\
        {
            b=(b>0 ? 1 : (b<0 ? -1 : 0));
        });
  }

  // 4. marching cubes�㷨�����������棺
  MatrixXd versResult_SDF, versResults_signs;
  MatrixXi trisResult_SDF, trisResults_signs;
  igl::marching_cubes(SDF, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResult_SDF, trisResult_SDF);
  igl::marching_cubes(signValues, gridCenters, gridCounts(0), gridCounts(1), gridCounts(2), 0, versResults_signs, trisResults_signs);


  // 5. ������ӻ���
  cout<<R"(Usage:
'1'  Show original mesh.
'2'  Show marching cubes contour of signed distance.
'3'  Show marching cubes contour of indicator function.
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(versResult_SDF, trisResult_SDF);
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
          viewer.data().clear();
          viewer.data().set_mesh(vers, tris);
          break;
        case '2':
          viewer.data().clear();
          viewer.data().set_mesh(versResult_SDF, trisResult_SDF);
          break;
        case '3':
          viewer.data().clear();
          viewer.data().set_mesh(versResults_signs, trisResults_signs);
          break;
      }
      viewer.data().set_face_based(true);
      return true;
    };

  viewer.launch();
}
