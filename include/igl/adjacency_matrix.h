#ifndef IGL_ADJACENCY_MATRIX_H
#define IGL_ADJACENCY_MATRIX_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  // Constructs the graph adjacency matrix  of a given mesh (vers,tris)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   tris  #tris by dim list of mesh simplices
  // Outputs: 
  //   A  max(tris)+1 by max(tris)+1 adjacency matrix, each row i corresponding to vers(i,:)
  //
  // Example:
  //   // Mesh in (vers,tris)
  //   Eigen::SparseMatrix<double> A;
  //   adjacency_matrix(tris,A);
  //   // sum each row 
  //   SparseVector<double> Asum;
  //   sum(A,1,Asum);
  //   // Convert row sums into diagonal of sparse matrix
  //   SparseMatrix<double> Adiag;
  //   diag(Asum,Adiag);
  //   // Build uniform laplacian
  //   SparseMatrix<double> U;
  //   U = A-Adiag;
  //
  // See also: edges, cotmatrix, diag
  template <typename DerivedF, typename T>
  IGL_INLINE void adjacency_matrix(
    const Eigen::MatrixBase<DerivedF> & tris, 
    Eigen::SparseMatrix<T>& A);


  // Constructs an vertex adjacency for a polygon mesh.
  //
  // Inputs:
  //   I  #I vectorized list of polygon corner indices into rows of some matrix vers
  //   C  #polygons+1 list of cumulative polygon sizes so that C(i+1)-C(i) =
  //     size of the ith polygon, and so I(C(i)) through I(C(i+1)-1) are the
  //     indices of the ith polygon
  // Outputs:
  //   A  max(I)+1 by max(I)+1 adjacency matrix, each row i corresponding to vers(i,:)
  //
  template <typename DerivedI, typename DerivedC, typename T>
  IGL_INLINE void adjacency_matrix(
    const Eigen::MatrixBase<DerivedI> & I,
    const Eigen::MatrixBase<DerivedC> & C,
    Eigen::SparseMatrix<T>& A);
}

#ifndef IGL_STATIC_LIBRARY
#  include "adjacency_matrix.cpp"
#endif

#endif
