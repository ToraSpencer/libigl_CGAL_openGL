#ifndef IGL_WINDING_NUMBER_H
#define IGL_WINDING_NUMBER_H
#include "igl_inline.h"
#include <Eigen/Core>


// Minimum number of iterms per openmp thread
#ifndef IGL_WINDING_NUMBER_OMP_MIN_VALUE
#  define IGL_WINDING_NUMBER_OMP_MIN_VALUE 1000
#endif
namespace igl
{
    // winding_number()——


    // 
  /*
       Computes the generalized winding number at each dim-dimensional query 
                    point in O with respect to the oriented  one-codimensional mesh (vers,tris). 

       This is equivalent to summing the subtended  signed angles/solid angles of each element in (vers,tris). 

       See, "Robust Inside-Outside Segmentation using Generalized Winding Numbers" [Jacobson et al. 2013].
  
       Inputs:
         vers           #vers by dim list of mesh vertex positions
         tris             #tris by dim list of mesh facets as indices into rows of vers. If dim==2,
                               then (vers,tris) describes a set of edges in the plane. If dim==3, then (vers,tris)
                               describes a triangle mesh/soup.
         O               #O by dim list of query points

       Output:
         wNums              #O by 1 list of winding numbers
  
       See also: igl::fast_winding_number
  */
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedO,
    typename DerivedW>
  IGL_INLINE void winding_number(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedO> & O,
    Eigen::PlainObjectBase<DerivedW> & wNums);


  //    计算单个顶点的缠绕数： 
  /* 
       Inputs:
        vers  n by dim list of vertex positions
        tris  #tris by dim list of triangle indices, minimum index is 0
        p  single origin position
       Outputs:
        w  winding number of this point
  */
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedp>
  IGL_INLINE typename DerivedV::Scalar winding_number(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<Derivedp> & p);
}

#ifndef IGL_STATIC_LIBRARY
#  include "winding_number.cpp"
#endif

#endif
