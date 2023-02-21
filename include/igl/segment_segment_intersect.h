#ifndef IGL_SEGMENT_SEGMENT_INTERSECT_H
#define IGL_SEGMENT_SEGMENT_INTERSECT_H


#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{

    // 线段求交——判断三维空间中两个线段是否相交；
    /*
     Determine whether two line segments A,B intersect
     A: p + t*r :  t \in [0,1]
     B: q + u*s :  u \in [0,1]

     Inputs:
           p  3-vector origin of segment A
           r  3-vector direction of segment A
           q  3-vector origin of segment B
           s  3-vector direction of segment B
          eps precision

     Outputs:
           t  scalar point of intersection along segment A, t \in [0,1]
           u  scalar point of intersection along segment B, u \in [0,1]
     Returns true if intersection
    */
    template<typename DerivedSource, typename DerivedDir>
    IGL_INLINE bool segment_segment_intersect(
            const Eigen::MatrixBase <DerivedSource> &p,
            const Eigen::MatrixBase <DerivedDir> &r,
            const Eigen::MatrixBase <DerivedSource> &q,
            const Eigen::MatrixBase <DerivedDir> &s,
            double &t,
            double &u,
            double eps = 1e-6
    );

}
#ifndef IGL_STATIC_LIBRARY
#   include "segment_segment_intersect.cpp"
#endif
#endif //IGL_SEGMENT_SEGMENT_INTERSECT_H
