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
         SDF   gridCountsX*gridCountsY*gridCountsZ list of values at each grid corner
             i.e. SDF(x + y*xres + z*xres*yres) for corner (x,y,z)
         gridCenters  gridCountsX*gridCountsY*gridCountsZ by 3 array of corresponding grid corner vertex locations
         gridCountsX  resolutions of the grid in x dimension
         gridCountsY  resolutions of the grid in y dimension
         gridCountsZ  resolutions of the grid in z dimension
         selectedSDF  the selectedSDF of the surface to reconstruct
       
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
    const Eigen::MatrixBase<DerivedS> & SDF,
    const Eigen::MatrixBase<DerivedGV> & gridCenters,
    const unsigned gridCountsX,
    const unsigned gridCountsY,
    const unsigned gridCountsZ,
    const typename DerivedS::Scalar selectedSDF,
    Eigen::PlainObjectBase<DerivedV> &versResult,
    Eigen::PlainObjectBase<DerivedF> &trisResult);


  template <
    typename DerivedS, 
    typename DerivedGV, 
    typename DerivedGI, 
    typename DerivedV, 
    typename DerivedF>
  IGL_INLINE void marching_cubes(
    const Eigen::MatrixBase<DerivedS> & SDF,
    const Eigen::MatrixBase<DerivedGV> & gridCenters,
    const Eigen::MatrixBase<DerivedGI> & GI,
    const typename DerivedS::Scalar selectedSDF,
    Eigen::PlainObjectBase<DerivedV> &versResult,
    Eigen::PlainObjectBase<DerivedF> &trisResult);
}

#ifndef IGL_STATIC_LIBRARY
#  include "marching_cubes.cpp"
#endif

#endif
