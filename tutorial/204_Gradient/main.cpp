#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>
#include "tutorial_shared_path.h"


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  MatrixXd vers;
  MatrixXi tris;
  VectorXd scalarField;	// һϵ�еı�������ֵ, scalarField == scalarField(vers);
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", vers, tris);
  igl::readDMAT(TUTORIAL_SHARED_PATH "/cheburashka-scalar.dmat", scalarField);		
	
  // 1. �����ݶȾ���G
  SparseMatrix<double> G;			// �ݶ�����
  igl::grad(vers, tris, G);

  // 2. ����scalarField���ݶ�
  MatrixXd GS = Map<const MatrixXd>((G*scalarField).eval().data(),tris.rows(),3);
  const VectorXd GS_mag = GS.rowwise().norm();			// Compute gradient magnitude

  // 3. �趨�ݶ�ָʾ�ߵ���󳤶�ȡ
  const double max_size = igl::avg_edge_length(vers, tris) / GS_mag.mean();

  // 4. ����ָʾ�ߡ���ÿ������Ƭ�����ϻ�һ��ָʾ�ߣ�����Ϊ�ݶȷ��� 
  MatrixXd baryCenters;
  igl::barycenter(vers, tris, baryCenters);
  const RowVector3d black(0,0,0);

  // 5. visualize:
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, tris);
  viewer.data().set_data(scalarField);
  viewer.data().add_edges(baryCenters, baryCenters+max_size*GS, black);
  viewer.data().show_lines = false;			// ����������
  viewer.launch();
}
