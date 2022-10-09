#ifndef IGL_COLLAPSE_EDGE_H
#define IGL_COLLAPSE_EDGE_H
#include "igl_inline.h"
#include "min_heap.h"
#include "decimate_callback_types.h"
#include <Eigen/Core>
#include <vector>
#include <set>


namespace igl
{
    // collpse_edge()——用于网格精简的边折叠算法；


    // 若uEdges, tris等数据中含有被折叠的边，将该行数据中的所有元素改写为IGL_COLLAPSE_EDGE_NULL；
  #define IGL_COLLAPSE_EDGE_NULL 0


    // 重载1.1——折叠一条边；
    /*
        要求输入网格是封闭的流形网格： 
    
    Assumes (vers,tris) is a closed manifold mesh (except for previously collapsed faces which should be set to:
          [IGL_COLLAPSE_EDGE_NULL IGL_COLLAPSE_EDGE_NULL IGL_COLLAPSE_EDGE_NULL].
     
     Collapses exactly two faces and exactly 3 uEdges from uEdges (edgeIdx and one side of each face gets collapsed to the other).

     This is implemented in a way that it can be repeatedly called until satisfaction and then the garbage in tris
          can be collected by removing NULL faces.

     Inputs:
       edgeIdx                      index into uEdges of edge to try to collapse. uEdges(edgeIdx,:) = [s d] or [d s] so that s<d, then d is collapsed to s.
       collapsedVer              dim list of vertex position where to place merged vertex

     Inputs/Outputs:
       vers                     by dim list of vertex positions, lesser index of uEdges(edgeIdx,:) will be set to midpoint of edge.
       tris                       by 3 list of face indices into vers.
       uEdges                  by 2 list of edge indices into vers.
       edgeUeInfo         tris*3 list of indices into uEdges, mapping each directed edge to unique unique edge in uEdges
       UeTrisInfo           #uEdges by 2 list of edge flaps, UeTrisInfo(edgeIdx,0)=f means edgeIdx=(i-->j) is the edge of
                                                tris(f,:) opposite the vth corner, where UeCornersInfo(edgeIdx,0)=v. 
                                                Similarly UeTrisInfo(edgeIdx,1) " edgeIdx=(j->i)
       UeCornersInfo    #uEdges by 2 list of edge flap corners (see above).
       e1                       index into uEdges of edge collpased on left
       e2                       index into uEdges of edge collpased on right
       f1                       index into tris of face collpased on left
       f2                       index into tris of face collpased on right

     Returns true if edge was collapsed
  */
  IGL_INLINE bool collapse_edge(
    const int edgeIdx,
    const Eigen::RowVectorXd & collapsedVer,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo,
    int & e1,
    int & e2,
    int & f1,
    int & f2);


  // 重载1——折叠一条边；
  IGL_INLINE bool collapse_edge(
    const int edgeIdx,
    const Eigen::RowVectorXd & collapsedVer,
    /*const*/ std::vector<int> & Nsv,
    const std::vector<int> & Nsf,
    /*const*/ std::vector<int> & Ndv,
    const std::vector<int> & Ndf,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo,
    int & e1,
    int & e2,
    int & f1,
    int & f2);


  // 重载2.1
  IGL_INLINE bool collapse_edge(
    const int edgeIdx,
    const Eigen::RowVectorXd & collapsedVer,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo);


  // 重载2.2
  /*
       Collapse least-cost edge from a priority queue and update queue 

       cost_and_placement()——function computing cost of collapsing an edge and 3d position where it should be placed:
       cost_and_placement( vers, tris, uEdges,
                        edgeUeInfo,UeTrisInfo,UeCornersInfo,
                        cost,
                        placement
                        );
           If the uEdges is collapsed
            then this function will be called on all uEdges of all faces previously incident on the endpoints of the collapsed edge.
  
       Inputs/Outputs:
         cost_and_placement         函数对象；
         Q                                       queue containing pairs of costs and edge indices and insertion "time"
         timeStamps                                     list of "time" of last time pushed into Q
         collapsedVers                   边折叠之后生成的顶点的坐标；list of stored placements
  */
  IGL_INLINE bool collapse_edge(
    const decimate_cost_and_placement_callback& cost_and_placement,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo,
    igl::min_heap< std::tuple<double,int,int>> & Q,
    Eigen::VectorXi & timeStamps,
    Eigen::MatrixXd & collapsedVers);


  // 重载2.3
  /*
   Inputs:
     pre_collapse  callback called with index of edge whose collapse is about
       to be attempted. This function should return whether to **proceed**
       with the collapse: returning true means "yes, try to collapse",
       returning false means "No, consider this edge 'uncollapsable', behave
       as if collapse_edge(edgeIdx) returned false.
     post_collapse  callback called with index of edge whose collapse was
       just attempted and a flag revealing whether this was successful.
  */
  IGL_INLINE bool collapse_edge(
    const decimate_cost_and_placement_callback & cost_and_placement,
    const decimate_pre_collapse_callback       & pre_collapse,
    const decimate_post_collapse_callback      & post_collapse,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo,
    igl::min_heap< std::tuple<double,int,int> > & Q,
    Eigen::VectorXi & timeStamps,
    Eigen::MatrixXd & collapsedVers);


  // 重载2
  /*
       Outputs:
         edgeIdx  index into uEdges of attempted collapsed edge. Set to -1 if Q is empty or
                contains only infinite cost uEdges.
         e1  index into uEdges of edge collpased on left.
         e2  index into uEdges of edge collpased on right.
         f1  index into tris of face collpased on left.
         f2  index into tris of face collpased on right.
  */
  IGL_INLINE bool collapse_edge(
    const decimate_cost_and_placement_callback & cost_and_placement,
    const decimate_pre_collapse_callback       & pre_collapse,
    const decimate_post_collapse_callback      & post_collapse,
    Eigen::MatrixXd & vers,
    Eigen::MatrixXi & tris,
    Eigen::MatrixXi & uEdges,
    Eigen::VectorXi & edgeUeInfo,
    Eigen::MatrixXi & UeTrisInfo,
    Eigen::MatrixXi & UeCornersInfo,
    igl::min_heap< std::tuple<double,int,int> > & Q,
    Eigen::VectorXi & timeStamps,
    Eigen::MatrixXd & collapsedVers,
    int & edgeIdx,
    int & e1,
    int & e2,
    int & f1,
    int & f2);
}

#ifndef IGL_STATIC_LIBRARY
#  include "collapse_edge.cpp"
#endif
#endif
