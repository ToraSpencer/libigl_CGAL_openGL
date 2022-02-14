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


// lambda���������¼�
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
        // ���¼�����������
        Eigen::SparseMatrix<double> mass;
        igl::massmatrix(newVers, tris, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

        // �����Է����� (mass - delta*L) * newVers = mass * newVers
        float delta = 0.001;
        const auto& S = (mass - delta * L);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
        assert(solver.info() == Eigen::Success);
        newVers = solver.solve(mass * newVers).eval();

#if 0
        // Compute centroid and subtract (also important for numerics)
        VectorXd dblA;                                     // ÿ������Ƭ�����������
        igl::doublearea(newVers, tris, dblA);
        double areaSum = 0.5 * dblA.sum();
        MatrixXd centers;
        igl::barycenter(newVers, tris, centers);
        RowVector3d centroid(0, 0, 0);
        for (int i = 0; i < centers.rows(); i++)
            centroid += 0.5 * dblA(i) / areaSum * centers.row(i);
        newVers.rowwise() -= centroid;

        // �����һ��
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

  // 1.a ֱ�ӹ���laplacian����Compute Laplace-Beltrami operator: 
  igl::cotmatrix(vers, tris, L);

  // 1.b �ֲ�����laplacian
  {
      SparseMatrix<double> Gradient, L2;

      igl::grad(vers, tris, Gradient);      // ��ɢ�ݶ�

      // Diagonal per-triangle "mass matrix"
      VectorXd dblA;
      igl::doublearea(vers, tris, dblA);             // ÿ������Ƭ���������

      // Place areas along diagonal  #dim times
      const auto& T = 1. * (dblA.replicate(3, 1) * 0.5).asDiagonal();
 
      L2 = -Gradient.transpose() * T * Gradient;         // discrete Dirichelet energy Hessian ��ɢ��������������������
      std::cout << "���ַ����õ���laplacian�Ĳ�ķ�����" << std::endl;
      cout << "(L2 - L).norm() == " << (L2 - L).norm() << endl;
  }

  // 2. ����ԭʼ�ķ�������ʹ��αɫ
  MatrixXd norms;
  igl::per_vertex_normals(vers, tris, norms);
  MatrixXd colors = norms.rowwise().normalized().array()*0.5+0.5;

  // 3. viewr����ʼ����
  newVers = vers;
  viewer.data().set_mesh(newVers, tris);
  viewer.data().set_colors(colors);
  viewer.callback_key_down = key_down;

  // 4. ����
  cout<<"Press [space] to smooth."<<endl;;
  cout<<"Press [r] to reset."<<endl;;
  return viewer.launch();
}
