#include "decimate_trivial_callbacks.h"

IGL_INLINE void igl::decimate_trivial_callbacks( decimate_pre_collapse_callback  & always_try,
    decimate_post_collapse_callback & never_care)
{
  always_try = [](
    const Eigen::MatrixXd &                             ,/*V*/
    const Eigen::MatrixXi &                             ,/*F*/
    const Eigen::MatrixXi &                             ,/*E*/
    const Eigen::VectorXi &                             ,/*EMAP*/
    const Eigen::MatrixXi &                             ,/*EF*/
    const Eigen::MatrixXi &                             ,/*EI*/
    const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
    const Eigen::VectorXi &                             ,/*EQ*/
    const Eigen::MatrixXd &                             ,/*C*/
    const int                                            /*e*/
    ) -> bool { return true;};


  never_care = [](
    const Eigen::MatrixXd &                             ,/*V*/
    const Eigen::MatrixXi &                             ,/*F*/
    const Eigen::MatrixXi &                             ,/*E*/
    const Eigen::VectorXi &                             ,/*EMAP*/
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
    )-> void { };
}
