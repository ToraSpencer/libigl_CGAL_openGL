#ifndef IGL_CIRCULATION_H
#define IGL_CIRCULATION_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
    // circulation()����Ѱ������ߵĶ˵����������Ƭ������1���򶥵㣻

    // ����1
    /*
       Return list of faces around the end point of an edge. 
       Assumes data-structures are built from an edge-manifold **closed** mesh.
  
       Inputs:
             edgeIdx               index into E of edge to circulate
             ccw           trueΪβ�˵㣬falseΪβ�˵㣻
             edgeUeInfo       #tris*3 list of indices into E, mapping each directed edge to unique unique edge in E
             EF             #E by 2 list of edge flaps, EF(edgeIdx,0)=f means edgeIdx=(i-->j) is the edge of
                                        tris(f,:) opposite the vth corner, where EI(edgeIdx,0)=v. Similarly EF(edgeIdx,1) "       edgeIdx=(j->i)
             EI              #E by 2 list of edge flap corners (see above).
        
        ��������Ƭ������
        Returns list of faces touched by circulation (in cyclically order).
     
       See also: edge_flaps
    */
  IGL_INLINE std::vector<int> circulation(
    const int edgeIdx,
    const bool ccw,
    const Eigen::VectorXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI);


  // ����1.1 Wrapper with VectorXi output.
  IGL_INLINE void circulation(
    const int edgeIdx,
    const bool ccw,
    const Eigen::VectorXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    Eigen::VectorXi & vN);


  // ����2��
  /*
       Outputs:
         nbrVersIdx            list of "next" vertex indices
         nbrTrisIdx            list of face indices
   */
  IGL_INLINE void circulation(
    const int edgeIdx,
    const bool ccw,
    const Eigen::MatrixXi & tris,
    const Eigen::VectorXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    std::vector<int> & nbrVersIdx,
    std::vector<int> & nbrTrisIdx);
}

#ifndef IGL_STATIC_LIBRARY
#  include "circulation.cpp"
#endif
#endif
