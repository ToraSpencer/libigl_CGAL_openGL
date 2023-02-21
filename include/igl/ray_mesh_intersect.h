#ifndef IGL_RAY_MESH_INTERSECT_H
#define IGL_RAY_MESH_INTERSECT_H
#include "igl_inline.h"
#include "Hit.h"
#include <Eigen/Core>
#include <vector>


namespace igl
{
    // 重载1：网格射线求交； 注：如果射线很多建议使用AABB.h
    /*
       Shoot a ray against a mesh (vers,tris) and collect all hits. 
       If you have many rays, consider using AABB.h
  
       Inputs:
         source  3-vector origin of ray
         dir  3-vector direction of ray
         vers  #vers by 3 list of mesh vertex positions
         tris  #tris by 3 list of mesh face indices into vers

       Outputs:
          hits  **sorted** list of hits
       Returns true if there were any hits (hits.size() > 0)
  
       See also: AABB.h
    */
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV, 
    typename DerivedF> 
  IGL_INLINE bool ray_mesh_intersect(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    std::vector<igl::Hit> & hits);


  // Outputs:
  //   hit  first hit, set only if it exists
  // Returns true if there was a hit
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV, 
    typename DerivedF> 
  IGL_INLINE bool ray_mesh_intersect(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    igl::Hit & hit);
}

#ifndef IGL_STATIC_LIBRARY
#  include "ray_mesh_intersect.cpp"
#endif
#endif
