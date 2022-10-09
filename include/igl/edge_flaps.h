#ifndef IGL_EDGE_FLAPS_H
#define IGL_EDGE_FLAPS_H
#include "igl_inline.h"
#include <Eigen/Core>


namespace igl
{

    // edge_flaps()
    /*
       Determine "edge flaps": two faces on either side of a unique edge (assumes edge-manifold mesh)
  
       Inputs:
         tris                            #tris by 3 list of face indices
         uEdges                     #uEdges by 2 list of edge indices into V.
         edgeUeInfo              #tris*3 list of indices into uEdges, mapping each directed edge to unique edge in uEdges

       Outputs:
         EF                         #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
                                               tris(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) " e=(j->i)
         EI                          #E by 2 list of edge flap corners (see above).
  
       See also: unique_edge_map
  
       TODO: 
            This seems to be a duplicate of edge_topology.h
            igl::edge_topology(V,tris,etEV,etFE,etEF);
            igl::edge_flaps(tris,efE,efEMAP,efEF,efEI);
            [~,I] = sort(efE,2)
            all( efE(sub2ind(size(efE),repmat(1:size(efE,1),2,1)',I)) == etEV )
            all( efEF(sub2ind(size(efE),repmat(1:size(efE,1),2,1)',I)) == etEF )
            all(efEMAP(sub2ind(size(tris),repmat(1:size(tris,1),3,1)',repmat([1 2 3],size(tris,1),1))) == etFE(:,[2 3 1]))
    */
  IGL_INLINE void edge_flaps(
    const Eigen::MatrixXi & tris,
    const Eigen::MatrixXi & uEdges,
    const Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI);


  // Only faces as input
  IGL_INLINE void edge_flaps(
    const Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI);
}


#ifndef IGL_STATIC_LIBRARY
#  include "edge_flaps.cpp"
#endif

#endif
