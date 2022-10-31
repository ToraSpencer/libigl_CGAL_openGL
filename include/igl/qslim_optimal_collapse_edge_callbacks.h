#ifndef IGL_QSLIM_OPTIMAL_COLLAPSE_EDGE_CALLBACKS_H
#define IGL_QSLIM_OPTIMAL_COLLAPSE_EDGE_CALLBACKS_H
#include "igl_inline.h"
#include "decimate_callback_types.h"
#include <Eigen/Core>
#include <functional>
#include <vector>
#include <tuple>
#include <set>


namespace igl
{
    // 生成qslim网格精简算法中需要使用的函数子—— cost_and_placement(), pre_collapse(), post_collapse();
    /*
       Prepare callbacks for decimating edges using the qslim optimal placement metric.
  
       Inputs:
         uEdges                       #uEdges by 2 list of working edges
         quadrics                       reference to list of working per vertex quadrics 
         v1                                  working variable to maintain end point of collapsed edge
         v2                                 working variable to maintain end point of collapsed edge

       Outputs
         cost_and_placement                 callback for evaluating cost of edge collapse and  determining placement of vertex (see collapse_edge)
         pre_collapse                               callback before edge collapse (see collapse_edge)
         post_collapse                              callback after edge collapse (see collapse_edge)
  */
  IGL_INLINE void qslim_optimal_collapse_edge_callbacks(
        Eigen::MatrixXi & uEdges,
        std::vector<std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double>>& quadrics,
        int & v1,
        int & v2,
        decimate_cost_and_placement_callback & cost_and_placement,
        decimate_pre_collapse_callback & pre_collapse,
        decimate_post_collapse_callback & post_collapse);
}

#ifndef IGL_STATIC_LIBRARY
#  include "qslim_optimal_collapse_edge_callbacks.cpp"
#endif
#endif
