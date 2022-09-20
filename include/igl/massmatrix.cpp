 #include "massmatrix.h"
#include "massmatrix_intrinsic.h"
#include "edge_lengths.h"
#include "normalize_row_sums.h"
#include "sparse.h"
#include "doublearea.h"
#include "repmat.h"
#include <Eigen/Geometry>
#include <iostream>



template <typename DerivedV, typename DerivedF, typename Scalar>
IGL_INLINE void igl::massmatrix(
  const Eigen::MatrixBase<DerivedV> & V, 
  const Eigen::MatrixBase<DerivedF> & F, 
  const MassMatrixType type,
  Eigen::SparseMatrix<Scalar>& M)
{
  using namespace Eigen;
  using namespace std;

  const int n = V.rows();
  const int m = F.rows();
  const int simplex_size = F.cols();

  MassMatrixType eff_type = type;
  // Use voronoi of for triangles by default, otherwise barycentric
  if(type == MASSMATRIX_TYPE_DEFAULT)
  {
    eff_type = (simplex_size == 3?MASSMATRIX_TYPE_VORONOI:MASSMATRIX_TYPE_BARYCENTRIC);
  }

  // Not yet supported
  assert(type!=MASSMATRIX_TYPE_FULL);

  if(simplex_size == 3)
  {
    // Triangles
    // edge lengths numbered same as opposite vertices
    Matrix<Scalar,Dynamic,3> l;
    igl::edge_lengths(V,F,l);
    return massmatrix_intrinsic(l,F,type,M);
  }else if(simplex_size == 4)
  {
    Matrix<typename DerivedF::Scalar,Dynamic,1> MI;
    Matrix<typename DerivedF::Scalar,Dynamic,1> MJ;
    Matrix<Scalar,Dynamic,1> MV;
    assert(V.cols() == 3);
    assert(eff_type == MASSMATRIX_TYPE_BARYCENTRIC);
    MI.resize(m*4,1); MJ.resize(m*4,1); MV.resize(m*4,1);
    MI.block(0*m,0,m,1) = F.col(0);
    MI.block(1*m,0,m,1) = F.col(1);
    MI.block(2*m,0,m,1) = F.col(2);
    MI.block(3*m,0,m,1) = F.col(3);
    MJ = MI;
    // loop over tets
    for(int i = 0;i<m;i++)
    {
      // http://en.wikipedia.org/wiki/Tetrahedron#Volume
      Matrix<Scalar,3,1> v0m3,v1m3,v2m3;
      v0m3.head(V.cols()) = V.row(F(i,0)) - V.row(F(i,3));
      v1m3.head(V.cols()) = V.row(F(i,1)) - V.row(F(i,3));
      v2m3.head(V.cols()) = V.row(F(i,2)) - V.row(F(i,3));
      Scalar v = fabs(v0m3.dot(v1m3.cross(v2m3)))/6.0;
      MV(i+0*m) = v/4.0;
      MV(i+1*m) = v/4.0;
      MV(i+2*m) = v/4.0;
      MV(i+3*m) = v/4.0;
    }
    sparse(MI,MJ,MV,n,n,M);
  }else
  {
    // Unsupported simplex size
    assert(false && "Unsupported simplex size");
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::massmatrix<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
// generated by autoexplicit.sh
template void igl::massmatrix<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
// generated by autoexplicit.sh
template void igl::massmatrix<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 4, 0, -1, 4>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 4, 0, -1, 4> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
// generated by autoexplicit.sh
template void igl::massmatrix<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
template void igl::massmatrix<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
template void igl::massmatrix<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MassMatrixType, Eigen::SparseMatrix<double, 0, int>&);
#endif
