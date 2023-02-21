#ifndef IGL_RAY_BOX_INTERSECT_H
#define IGL_RAY_BOX_INTERSECT_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>


namespace igl
{
    // 射线-盒求交——判断射线和一个轴向和是否相交；
    /*
       Determine whether a ray origin+t*dir and box intersect within the ray's parameterized range (t0,t1)
  
       Inputs:
             source  3-vector origin of ray
             dir  3-vector direction of ray
             box  axis aligned box
             t0  hit only if hit.t less than t0
             t1  hit only if hit.t greater than t1

       Outputs:
             tmin  minimum of interval of overlap within [t0,t1]
             tmax  maximum of interval of overlap within [t0,t1]

       Returns true if hit
    */
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename Scalar>
  IGL_INLINE bool ray_box_intersect(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::AlignedBox<Scalar,3> & box,
    const Scalar & t0,
    const Scalar & t1,
    Scalar & tmin,
    Scalar & tmax);
}
#ifndef IGL_STATIC_LIBRARY
#  include "ray_box_intersect.cpp"
#endif
#endif
