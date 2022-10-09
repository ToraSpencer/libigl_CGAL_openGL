 #ifndef IGL_DECIMATE_H
#define IGL_DECIMATE_H
#include "igl_inline.h"
#include "decimate_callback_types.h"
#include <Eigen/Core>

namespace igl
{
    // 使用边折叠算法精简网格：

    // 重载1.1
      /*
           Assumes (vers,tris) is a manifold mesh (possibly with boundary) Collapses edges
           until desired number of faces is achieved. This uses default edge cost and
           merged vertex placement functions {edge length, edge midpoint}.
  
           Inputs:
             vers                       #vers by dim list of vertex positions
             tris                        #tris by 3 list of face indices into vers.
             maxTrisCount                  desired number of output faces

           Outputs:
             versOut                            #versOut by dim list of output vertex posistions (can be same ref as vers)
             trisOut                            #trisOut by 3 list of output face indices into versOut (can be same ref as trisOut)
             newOldTrisInfo                 #trisOut list of indices into tris of birth face
             newOldVersInfo                 #versOut list of indices into vers of birth vertices

           Returns true if m was reached (otherwise #trisOut > m)
        */
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const size_t maxTrisCount,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo,
    Eigen::VectorXi & newOldVersInfo);


  // 重载1.1。1
  /*
       Inputs:
         vers  #vers by dim list of vertex positions
         tris  #tris by 3 list of face indices into vers.
         maxTrisCount  desired number of output faces

       Outputs:
         versOut  #versOut by dim list of output vertex posistions (can be same ref as vers)
         trisOut  #trisOut by 3 list of output face indices into versOut (can be same ref as trisOut)
         newOldTrisInfo  #trisOut list of indices into tris of birth face
       Returns true if m was reached (otherwise #trisOut > m)
  */
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const size_t maxTrisCount,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo);


  // 重载1.2——使用用户指定的边cost计算方法，以及循环终止条件；
  /*
       Assumes a **closed** manifold mesh. 
       See igl::connect_boundary_to_infinity  and igl::decimate in decimate.cpp

       is handling meshes with boundary by connecting all boundary edges with dummy facets 
            to infinity and modifying the stopping criteria.
  

       Inputs:
         cost_and_placement              function computing cost of collapsing an edge and 3d
                                                               position where it should be placed: cost_and_placement(vers,tris,edges,edgeUeInfo,EF,EI,cost,placement);
         
         stopping_condition                function returning whether to stop collapsing edges
                                                               based on current state. Guaranteed to be called after _successfully_
                                                               collapsing edge e removing edges (e,e1,e2) and faces (f1,f2):
                                                               bool should_stop =  stopping_condition(vers,tris,edges,edgeUeInfo,EF,EI,Q,Qit,C,e,e1,e2,f1,f2);
    */
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const decimate_cost_and_placement_callback & cost_and_placement,
    const decimate_stopping_condition_callback & stopping_condition,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo,
    Eigen::VectorXi & newOldVersInfo);


  // 重载1.3
  /*
       Inputs:
         pre_collapse               callback called with index of edge whose collapse is about  to be attempted (see collapse_edge)

         post_collapse              callback called with index of edge whose collapse was
                                                   just attempted and a flag revealing whether this was successful (see collapse_edge)
  */
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const decimate_cost_and_placement_callback & cost_and_placement,
    const decimate_stopping_condition_callback & stopping_condition,
    const decimate_pre_collapse_callback       & pre_collapse,
    const decimate_post_collapse_callback      & post_collapse,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo,
    Eigen::VectorXi & newOldVersInfo);


  // 重载1.
  /*
       Inputs:
         edgeUeInfo         #tris*3 list of indices into edges, mapping each directed edge to unique  unique edge in edges

         EF                         #edges by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
                                            tris(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "  e=(j->i)

         EI                         #edges by 2 list of edge flap corners (see above).
  */
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const decimate_cost_and_placement_callback & cost_and_placement,
    const decimate_stopping_condition_callback & stopping_condition,
    const decimate_pre_collapse_callback       & pre_collapse,
    const decimate_post_collapse_callback      & post_collapse,
    const Eigen::MatrixXi & edges,
    const Eigen::VectorXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo,
    Eigen::VectorXi & newOldVersInfo);

}

#ifndef IGL_STATIC_LIBRARY
#  include "decimate.cpp"
#endif
#endif


