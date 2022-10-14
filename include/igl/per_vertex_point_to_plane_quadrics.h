#ifndef IGL_PER_VERTEX_POINT_TO_PLANE_QUADRICS_H
#define IGL_PER_VERTEX_POINT_TO_PLANE_QUADRICS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <tuple>


namespace igl
{

    // per_vertex_point_to_plane_quadrics()¡ª¡ª

    /*
       Compute quadrics per vertex of a "closed" triangle mesh (vers,tris). 

       Rather than follow the qslim paper, this implements the lesser-known _follow up_
            "Simplifying Surfaces with Color and Texture using Quadric Error Metrics".

       This allows vers to be n-dimensional (where the extra coordiantes store texture UVs, color RGBs, etc.
  
       Inputs:
         vers                   #vers by n list of vertex positions. Assumes that vertices with
                                           infinite coordinates are "points at infinity" being used to close up
                                           boundary edges with faces. This allows special subspace quadrice for
                                           boundary edges: There should never be more than one "point at infinity" in a single triangle.
         tris                   #tris by 3 list of triangle indices into vers
         edgeUeInfo       #tris*3 list of indices into E, mapping each directed edge to unique unique edge in E
         EF                     #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
                                         tris(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "  e=(j->i)
         EI                     #E by 2 list of edge flap corners (see above).

       Outputs:
         quadrics   #vers list of quadrics, where a quadric is a tuple {A,b,c} such that the quadratic energy of 
                                    moving this vertex to position x is given by x'Ax - 2b + c
  */
  IGL_INLINE void per_vertex_point_to_plane_quadrics(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const Eigen::MatrixXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    std::vector<
      std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> > & quadrics);
}
#ifndef IGL_STATIC_LIBRARY
#  include "per_vertex_point_to_plane_quadrics.cpp"
#endif
#endif
