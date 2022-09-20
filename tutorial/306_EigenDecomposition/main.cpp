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

// ��������LB���ӵ������ֽ�
/*
     LB(vers) == M.inverse() * L * vers;
     �� LB*x == lambda * x;     �ȼ��� L*x == lambda* M.inverse() * x;
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
  SparseMatrix<double> L, M;            // LB���ӣ� ������������
  igl::opengl::glfw::Viewer viewer;

  if(!read_triangle_mesh(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off", vers, tris))
        cout<<"failed to load mesh"<<endl;

  // 1. 
  is2D = vers.col(2).minCoeff()==vers.col(2).maxCoeff();
  scale = (vers.colwise().maxCoeff()-vers.colwise().minCoeff()).norm();

  // 2. ����LB���ӣ�������������
  cotmatrix(vers, tris, L);
  L = (-L).eval();
  massmatrix(vers, tris, MASSMATRIX_TYPE_DEFAULT, M);

  // 3. ��������ֵ����������
  const size_t k = 5;
  if(!eigs(L, M, k+1, EIGS_TYPE_SM, eigVecs, eigValues))
        cout<<"failed."<<endl;
  eigVecs = ((eigVecs.array()-eigVecs.minCoeff())/(eigVecs.maxCoeff()-eigVecs.minCoeff())).eval();  // Normalize


  // lambda���������¼���Ӧ��
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int)->bool
  {
    switch(key)
    {
      default:
        return false;
      case ' ':
      {
          // ���¿ո��ʱ��ʾ��ͬeigValue-eigVec�Ե����ݣ�
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
