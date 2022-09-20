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
  VectorXd scalarField;	// 一系列的标量函数值, scalarField == scalarField(vers);
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", vers, tris);
  igl::readDMAT(TUTORIAL_SHARED_PATH "/cheburashka-scalar.dmat", scalarField);		
	
  // 1. 计算梯度矩阵G
  SparseMatrix<double> G;			// 梯度算子
  igl::grad(vers, tris, G);

  // 2. 计算scalarField的梯度
  MatrixXd GS = Map<const MatrixXd>((G*scalarField).eval().data(),tris.rows(),3);
  const VectorXd GS_mag = GS.rowwise().norm();			// Compute gradient magnitude

  // 3. 设定梯度指示线的最大长度取
  const double max_size = igl::avg_edge_length(vers, tris) / GS_mag.mean();

  // 4. 计算指示线——每个三角片重心上画一根指示线，方向为梯度方向。 
  MatrixXd baryCenters;
  igl::barycenter(vers, tris, baryCenters);
  const RowVector3d black(0,0,0);

  // 5. visualize:
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(vers, tris);
  viewer.data().set_data(scalarField);
  viewer.data().add_edges(baryCenters, baryCenters+max_size*GS, black);
  viewer.data().show_lines = false;			// 隐藏网格线
  viewer.launch();
}
