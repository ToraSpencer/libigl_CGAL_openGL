#ifndef IGL_EDGE_FLAPS_H
#define IGL_EDGE_FLAPS_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
    // edge_flaps()——计算每条无向边关联的三角片，默认输入的是流形网格（最多只关联两个三角片，边缘边关联一个三角片，另一个索引写-1）
    /*
      定义角(corner)的索引，三角片第i个顶点所在的角的索引为i；
    
       Determine "edge flaps": two faces on either side of a unique edge (assumes edge-manifold mesh)
  
       Inputs:
         tris                            #tris by 3 list of face indices
         uEdges                     #uEdges by 2 list of edge indices into vers.
         edgeUeInfo              #tris*3 list of indices into uEdges, mapping each directed edge to unique edge in uEdges

       Outputs:
         UeTrisInfo                第i行的元素：索引为i的无向边关联的三角片的索引；
                                            list of edge flaps, UeTrisInfo(e,0) = f means e=(i-->j) is the edge of tris(f,:) opposite the vth corner, 
                                                 where UeCornersInfo(e,0)=v. Similarly UeTrisInfo(e,1) " e=(j->i)
         UeCornersInfo          边相对的两个角的索引，不存在的话写-1
                                          list of edge flap corners (see above).
  
       See also: unique_edge_map
  
       TODO: 
            This seems to be a duplicate of edge_topology.h
            igl::edge_topology(vers,tris,etEV,etFE,etEF);
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
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo);


  // Only faces as input
  IGL_INLINE void edge_flaps(
    const Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo);
}


#ifndef IGL_STATIC_LIBRARY
#  include "edge_flaps.cpp"
#endif

#endif
