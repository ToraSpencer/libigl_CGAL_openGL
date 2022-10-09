 #ifndef IGL_MASSMATRIX_H
#define IGL_MASSMATRIX_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>


// 计算质量矩阵
namespace igl 
{

  enum MassMatrixType
  {
    MASSMATRIX_TYPE_BARYCENTRIC = 0,
    MASSMATRIX_TYPE_VORONOI = 1,
    MASSMATRIX_TYPE_FULL = 2,
    MASSMATRIX_TYPE_DEFAULT = 3,
    NUM_MASSMATRIX_TYPE = 4
  };


  /*
   Constructs the mass (area) matrix for a given mesh (vers,tris).
  
       Templates:
         DerivedV  derived type of eigen matrix for vers (e.g. derived from
           MatrixXd)
         DerivedF  derived type of eigen matrix for tris (e.g. derived from
           MatrixXi)
         Scalar  scalar type for eigen sparse matrix (e.g. double)

       Inputs:
         vers  #vers by dim list of mesh vertex positions
         tris  #tris by simplex_size list of mesh elements (triangles or tetrahedra)
         type  one of the following ints:
               MASSMATRIX_TYPE_BARYCENTRIC  barycentric——M(i, i)是顶点i的重心区域的面积；
               MASSMATRIX_TYPE_VORONOI voronoi-hybrid {default}——M(i, i)是顶点i的voronoi区域的面积；
               MASSMATRIX_TYPE_FULL full {not implemented}

       Outputs: 
         M  #vers by #vers mass matrix

  */
  template <typename DerivedV, typename DerivedF, typename Scalar>
  IGL_INLINE void massmatrix(
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    const MassMatrixType type,
    Eigen::SparseMatrix<Scalar>& M);
}

#ifndef IGL_STATIC_LIBRARY
#  include "massmatrix.cpp"
#endif

#endif

