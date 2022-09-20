 
#ifndef IGL_COTMATRIX_H
#define IGL_COTMATRIX_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

// 生成余切权的laplacian (同时也是刚度矩阵(stiffness matrix))

// History:
/*  
    Used const references rather than copying the entire mesh
        Alec 9 October 2011
      removed cotan (uniform weights) optional parameter it was building a buggy
        half of the uniform laplacian, please see adjacency_matrix instead 
     Alec 9 October 2011
*/   

namespace igl 
{
  /*
      Constructs the cotangent stiffness matrix(discrete laplacian) for a given mesh (vers,tris).
  
       Templates:
         DerivedV  derived type of eigen matrix for vers (e.g. derived from
           MatrixXd)
         DerivedF  derived type of eigen matrix for tris (e.g. derived from
           MatrixXi)
         Scalar  scalar type for eigen sparse matrix (e.g. double)

       Inputs:
         vers  #vers by dim list of mesh vertex positions
         tris  #tris by simplex_size list of mesh elements (triangles or tetrahedra)

       Outputs: 
         L  #vers by #vers cotangent matrix, each row i corresponding to vers(i,:)
  
       See also: adjacency_matrix
  
       Note: This Laplacian uses the convention that diagonal entries are
       **minus** the sum of off-diagonal entries. The diagonal entries are
       therefore in general negative and the matrix is **negative** semi-definite
       (immediately, -L is **positive** semi-definite)
  */
  template <typename DerivedV, typename DerivedF, typename Scalar>
  IGL_INLINE void cotmatrix(
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    Eigen::SparseMatrix<Scalar>& L);


  /*
   Cotangent Laplacian(and mass matrix) for polygon meshes according to
       "Polygon Laplacian Made Simple" [Bunge et al. 2020]
  
       Inputs:
         vers  #vers by 3 list of mesh vertex positions
         I  #I vectorized list of polygon corner indices into rows of some matrix vers
         C  #polygons+1 list of cumulative polygon sizes so that C(i+1)-C(i) = size of
           the ith polygon, and so I(C(i)) through I(C(i+1)-1) are the indices of
           the ith polygon
   
       Outputs:
         L  #vers by #vers polygon Laplacian made simple matrix
         M  #vers by #vers mass matrix
         P  #vers+#polygons by #vers prolongation operator
  */

  template <
    typename DerivedV, 
    typename DerivedI, 
    typename DerivedC, 
    typename Scalar>
  IGL_INLINE void cotmatrix(
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedI> & I, 
    const Eigen::MatrixBase<DerivedC> & C, 
    Eigen::SparseMatrix<Scalar>& L,
    Eigen::SparseMatrix<Scalar>& M,
    Eigen::SparseMatrix<Scalar>& P);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cotmatrix.cpp"
#endif

#endif
