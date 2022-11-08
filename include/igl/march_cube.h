#ifndef IGL_MARCH_CUBE_H
#define IGL_MARCH_CUBE_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <unordered_map>


namespace igl
{
    // march_cube()——生成marching cubes算法中的单个方块；
  /*
        Process a single cube of a marching cubes grid.
  
           Inputs:
             gridCenters          #gridCenters by 3 list of grid vertex positions
             cS                     list of 8 scalar field values at grid corners
             cI                         list of 8 indices of corners into rows of gridCenters
             isovalue        level-set value being extracted (often 0)
             vers          #vers by 3 current list of output mesh vertex positions
             n                  current number of mesh vertices (i.e., occupied rows in vers)
             tris               #tris by 3 current list of output mesh triangle indices into rows of vers
             m                  current number of mesh triangles (i.e., occupied rows in tris)
             E2V        current edge (GV_i,GV_j) to vertex (V_k) map

           Side-effects: vers,n,tris,m,E2V are updated to contain new vertices and faces of any constructed mesh elements
  */
  template <
    typename DerivedGV,
    typename Scalar,
    typename Index,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE void march_cube(
    const DerivedGV & gridCenters,
    const Eigen::Matrix<Scalar,8,1> & cS,
    const Eigen::Matrix<Index,8,1> & cI,
    const Scalar & isovalue,
    Eigen::PlainObjectBase<DerivedV> &vers,
    Index & n,
    Eigen::PlainObjectBase<DerivedF> &tris,
    Index & m,
    std::unordered_map<int64_t,int> & E2V);
}

#ifndef IGL_STATIC_LIBRARY
#  include "march_cube.cpp"
#endif

#endif
