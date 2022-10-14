#ifndef IGL_SHORTEST_EDGE_AND_MIDPOINT_H
#define IGL_SHORTEST_EDGE_AND_MIDPOINT_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{
    // shortest_edge_and_midpoint()�������۵��㷨�м���Ĭ�ϵ�cost���Լ��۵���Ķ��㣻
    /*
            ʹ�ñ߳���Ϊ���۵���cost���۵���Ķ���ȡ�ߵ��е㣻
           Cost and placement function compatible with igl::decimate. 
           The cost of collapsing an edge is its length (prefer to collapse short edges) 
           and the placement strategy for the new vertex is the midpoint of the collapsed edge.
    
           Inputs:
             edgeIdx        index into uEdges of edge to be considered for collapse
             vers               #vers by dim list of vertex positions
             tris               #tris by 3 list of faces (ignored)
             uEdges         #uEdges by 2 list of edge indices into vers
             edgeUeInfo           #tris*3 list of half-edges indices into uEdges (ignored)
             EF                 #uEdges by 2 list of edge-face flaps into tris (ignored)
             EI                 #uEdges by 2 list of edge-face opposite corners (ignored)

           Outputs:
             cost                           set to edge length
             edgeCenter              placed point set to edge midpoint
    */
  IGL_INLINE void shortest_edge_and_midpoint(
    const int edge,
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & /*tris*/,
    const Eigen::MatrixXi & uEdges,
    const Eigen::VectorXi & /*edgeUeInfo*/,
    const Eigen::MatrixXi & /*EF*/,
    const Eigen::MatrixXi & /*EI*/,
    double & cost,
    Eigen::RowVectorXd & edgeCenter);
}

#ifndef IGL_STATIC_LIBRARY
#  include "shortest_edge_and_midpoint.cpp"
#endif
#endif


