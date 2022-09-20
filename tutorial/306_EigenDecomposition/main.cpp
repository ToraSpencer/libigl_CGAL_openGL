#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>
#include "tutorial_shared_path.h"

// 求网格上LB算子的特征分解
/*
     LB(vers) == M.inverse() * L * vers;
     → LB*x == lambda * x;     等价于 L*x == lambda* M.inverse() * x;
*/


Eigen::MatrixXd vers, eigVecs;
Eigen::MatrixXi tris;
int index = 0;
double scale = 1;
bool is2D = 0;


int main(int argc,  char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  VectorXd eigValues;
  SparseMatrix<double> L, M;            // LB算子， 顶点质量矩阵
  igl::opengl::glfw::Viewer viewer;

  if(!read_triangle_mesh(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off", vers, tris))
        cout<<"failed to load mesh"<<endl;

  // 1. 
  is2D = vers.col(2).minCoeff()==vers.col(2).maxCoeff();
  scale = (vers.colwise().maxCoeff()-vers.colwise().minCoeff()).norm();

  // 2. 计算LB算子，顶点质量矩阵：
  cotmatrix(vers, tris, L);
  L = (-L).eval();
  massmatrix(vers, tris, MASSMATRIX_TYPE_DEFAULT, M);

  // 3. 计算特征值特征向量；
  const size_t k = 5;
  if(!eigs(L, M, k+1, EIGS_TYPE_SM, eigVecs, eigValues))
        cout<<"failed."<<endl;
  eigVecs = ((eigVecs.array()-eigVecs.minCoeff())/(eigVecs.maxCoeff()-eigVecs.minCoeff())).eval();  // Normalize


  // lambda——键盘事件相应；
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int)->bool
  {
    switch(key)
    {
      default:
        return false;
      case ' ':
      {
          // 按下空格键时显示不同eigValue-eigVec对的数据；
        eigVecs = eigVecs.rightCols(k).eval();
        VectorXd Z = scale*0.5*eigVecs.col(index);            // Rescale eigen vectors for visualization

        if(is2D)
        {
          vers.col(2) = Z;
          viewer.data().set_mesh(vers, tris);
          viewer.data().compute_normals();
        }

        viewer.data().set_data(eigVecs.col(index).eval());
        index = (index+1)%eigVecs.cols();
        return true;
      }
    }
  };

  // 4. visualize:
  viewer.data().set_mesh(vers, tris);
  viewer.callback_key_down(viewer, ' ', 0);
  viewer.data().show_lines = false;
  std::cout<< R"( [space] Cycle through eigen modes)";
  viewer.launch();
}
