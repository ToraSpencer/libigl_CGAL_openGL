#ifndef IGL_DECIMATE_CALLBACK_TYPES_H
#define IGL_DECIMATE_CALLBACK_TYPES_H
#include <Eigen/Core>
#include "min_heap.h"

// 定义了各种函数包装器std::function的别名；
namespace igl
{
  // Function handles used to customize the `igl::decimate` command.
  using decimate_cost_and_placement_callback = 
    std::function<void(
      const int                                           ,/*e*/
      const Eigen::MatrixXd &                             ,/*vers*/
      const Eigen::MatrixXi &                             ,/*tris*/
      const Eigen::MatrixXi &                             ,/*E*/
      const Eigen::VectorXi &                             ,/*edgeUeInfo*/
      const Eigen::MatrixXi &                             ,/*EF*/
      const Eigen::MatrixXi &                             ,/*EI*/
      double &                                            ,/*cost*/
      Eigen::RowVectorXd &                                 /*p*/
      )>;

  using decimate_stopping_condition_callback = 
    std::function<bool(
      const Eigen::MatrixXd &                             ,/*vers*/
      const Eigen::MatrixXi &                             ,/*tris*/
      const Eigen::MatrixXi &                             ,/*E*/
      const Eigen::VectorXi &                             ,/*edgeUeInfo*/
      const Eigen::MatrixXi &                             ,/*EF*/
      const Eigen::MatrixXi &                             ,/*EI*/
      const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
      const Eigen::VectorXi &                             ,/*EQ*/
      const Eigen::MatrixXd &                             ,/*C*/
      const int                                           ,/*e*/
      const int                                           ,/*e1*/
      const int                                           ,/*e2*/
      const int                                           ,/*f1*/
      const int                                            /*f2*/
      )>;

  using decimate_pre_collapse_callback = 
    std::function<bool(
      const Eigen::MatrixXd &                             ,/*vers*/
      const Eigen::MatrixXi &                             ,/*tris*/
      const Eigen::MatrixXi &                             ,/*E*/
      const Eigen::VectorXi &                             ,/*edgeUeInfo*/
      const Eigen::MatrixXi &                             ,/*EF*/
      const Eigen::MatrixXi &                             ,/*EI*/
      const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
      const Eigen::VectorXi &                             ,/*EQ*/
      const Eigen::MatrixXd &                             ,/*C*/
      const int                                            /*e*/
      )>;

  using decimate_post_collapse_callback = 
    std::function<void(
      const Eigen::MatrixXd &                             ,/*vers*/
      const Eigen::MatrixXi &                             ,/*tris*/
      const Eigen::MatrixXi &                             ,/*E*/
      const Eigen::VectorXi &                             ,/*edgeUeInfo*/
      const Eigen::MatrixXi &                             ,/*EF*/
      const Eigen::MatrixXi &                             ,/*EI*/
      const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
      const Eigen::VectorXi &                             ,/*EQ*/
      const Eigen::MatrixXd &                             ,/*C*/
      const int                                           ,/*e*/
      const int                                           ,/*e1*/
      const int                                           ,/*e2*/
      const int                                           ,/*f1*/
      const int                                           ,/*f2*/
      const bool                                           /*collapsed*/
      )>;
}
#endif
