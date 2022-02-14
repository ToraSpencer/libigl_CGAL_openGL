#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/repdiag.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd vers, newVers;
Eigen::MatrixXi tris;
Eigen::SparseMatrix<double> L;
igl::opengl::glfw::Viewer viewer;


// lambda——键盘事件
const auto& key_down = [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
{
    switch (key)
    {
    case 'r':

    case 'R':
        newVers = vers;
        break;

    case ' ':
    {
        // 重新计算质量矩阵
        Eigen::SparseMatrix<double> mass;
        igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

        // 解线性方程组 (mass - delta*L) * newVers = mass * newVers
        float delta = 0.001;
        const auto& S = (mass - delta * L);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
        assert(solver.info() == Eigen::Success);
        newVers = solver.solve(mass * newVers).eval();

#if 0
        // Compute centroid and subtract (also important for numerics)
        VectorXd dblA;                                     // 每个三角片面积的两倍；
        igl::doublearea(newVers, tris, dblA);
        double areaSum = 0.5 * dblA.sum();
        MatrixXd centers;
        igl::barycenter(newVers, tris, centers);
        RowVector3d centroid(0, 0, 0);
        for (int i = 0; i < centers.rows(); i++)
            centroid += 0.5 * dblA(i) / areaSum * centers.row(i);
        newVers.rowwise() -= centroid;

        // 面积归一化
        newVers.array() /= sqrt(areaSum);
#endif

        break;
    }

    default:
        return false;
    }

    viewer.data().set_vertices(newVers);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(newVers, tris);
    return true;
};


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", vers, tris);

  // 1.a 直接构造laplacian——Compute Laplace-Beltrami operator: 
  igl::cotmatrix(vers, tris, L);

  // 1.b 分步构造laplacian
  {
      SparseMatrix<double> Gradient, L2;

      igl::grad(vers, tris, Gradient);      // 离散梯度

      // Diagonal per-triangle "mass matrix"
      VectorXd dblA;
      igl::doublearea(vers, tris, dblA);             // 每个三角片面积的两倍

      // Place areas along diagonal  #dim times
      const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();
 
      L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian 离散狄利克雷能量海塞矩阵
      std::cout << "两种方法得到的laplacian的差的范数：" << std::endl;
      cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
  }

  // 2. 根据原始的法向量，使用伪色
  MatrixXd norms;
  igl::per_vertex_normals(vers, tris, norms);
  MatrixXd colors = norms.rowwise().normalized().array()*0.5+0.5;

  // 3. viewr填充初始数据
  newVers = vers;
  viewer.data().set_mesh(newVers, tris);
  viewer.data().set_colors(colors);
  viewer.callback_key_down = key_down;

  // 4. 运行
  cout<<"Press [space] to smooth."<<endl;;
  cout<<"Press [r] to reset."<<endl;;
  return viewer.launch();
}
