#ifndef IGL_MARCHING_CUBES_H
#define IGL_MARCHING_CUBES_H
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl
{
  /*
       performs marching cubes reconstruction on a grid defined by values, and points, 
            and generates a mesh defined by vertices and faces
  
       Input:
         scalarField   nx*ny*nz list of values at each grid corner
             i.e. scalarField(x + y*xres + z*xres*yres) for corner (x,y,z)
         gridCenters  nx*ny*nz by 3 array of corresponding grid corner vertex locations
         nx  resolutions of the grid in x dimension
         ny  resolutions of the grid in y dimension
         nz  resolutions of the grid in z dimension
         isovalue  the isovalue of the surface to reconstruct
       
       Output:
         versResult  #versResult by 3 list of mesh vertex positions
         trisResult  #trisResult by 3 list of mesh triangle indices into rows of versResult
  */
  template <
    typename DerivedS, 
    typename DerivedGV, 
    typename DerivedV, 
    typename DerivedF>
  IGL_INLINE void marching_cubes(
    const Eigen::MatrixBase<DerivedS> & scalarField,
    const Eigen::MatrixBase<DerivedGV> & gridCenters,
    const unsigned nx,
    const unsigned ny,
    const unsigned nz,
    const typename DerivedS::Scalar isovalue,
    Eigen::PlainObjectBase<DerivedV> &versResult,
    Eigen::PlainObjectBase<DerivedF> &trisResult);


  template <
    typename DerivedS, 
    typename DerivedGV, 
    typename DerivedGI, 
    typename DerivedV, 
    typename DerivedF>
  IGL_INLINE void marching_cubes(
    const Eigen::MatrixBase<DerivedS> & scalarField,
    const Eigen::MatrixBase<DerivedGV> & gridCenters,
    const Eigen::MatrixBase<DerivedGI> & GI,
    const typename DerivedS::Scalar isovalue,
    Eigen::PlainObjectBase<DerivedV> &versResult,
    Eigen::PlainObjectBase<DerivedF> &trisResult);
}

#ifndef IGL_STATIC_LIBRARY
#  include "marching_cubes.cpp"
#endif

#endif
