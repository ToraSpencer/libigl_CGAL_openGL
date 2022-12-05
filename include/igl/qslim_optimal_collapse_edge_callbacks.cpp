#include "qslim_optimal_collapse_edge_callbacks.h"
#include "quadric_binary_plus_operator.h"
#include <Eigen/LU>


// 生成qslim网格精简算法中需要使用的函数子—— cost_and_placement(), pre_collapse(), post_collapse();
IGL_INLINE void igl::qslim_optimal_collapse_edge_callbacks(
  Eigen::MatrixXi & uEdges,
  std::vector<std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> > &quadrics,
  int & v1,
  int & v2,
  decimate_cost_and_placement_callback & cost_and_placement,
  decimate_pre_collapse_callback       & pre_collapse,
  decimate_post_collapse_callback      & post_collapse)
{
  typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;

  // lambda——cost_and_placement()——qslim算法中计算每条边折叠的cost值，以及折叠后顶点的位置：
  cost_and_placement = [&quadrics,&v1,&v2](
    const int ueIdx,
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & /*tris*/,
    const Eigen::MatrixXi & uEdges,
    const Eigen::VectorXi & /*edgeUeInfo*/,
    const Eigen::MatrixXi & /*EF*/,
    const Eigen::MatrixXi & /*EI*/,
    double & cost,
    Eigen::RowVectorXd & place)
  {
    // Combined quadric
    Quadric quadric_p;
    quadric_p = quadrics[uEdges(ueIdx,0)] + quadrics[uEdges(ueIdx,1)];

    // Quadric: place'Ap + 2b'place + c,  optimal point: Ap = -b, or rather because we have row vectors: pA=-b
    const auto & A = std::get<0>(quadric_p);
    const auto & b = std::get<1>(quadric_p);
    const auto & c = std::get<2>(quadric_p);
    place = -b*A.inverse();
    cost = place.dot(place*A) + 2*place.dot(b) + c;

    // Force infs and nans to infinity
    if(std::isinf(cost) || cost!=cost)
    {
      cost = std::numeric_limits<double>::infinity();
      // Prevent NaNs. Actually NaNs might be useful for debugging.
      place.setConstant(0);
    }
  };


  // lambda——pre_collapse()
  pre_collapse = [&v1,&v2](
    const Eigen::MatrixXd &                             ,/*vers*/
    const Eigen::MatrixXi &                             ,/*tris*/
    const Eigen::MatrixXi & uEdges                           ,
    const Eigen::VectorXi &                             ,/*edgeUeInfo*/
    const Eigen::MatrixXi &                             ,/*EF*/
    const Eigen::MatrixXi &                             ,/*EI*/
    const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
    const Eigen::VectorXi &                             ,/*EQ*/
    const Eigen::MatrixXd &                             ,/*C*/
    const int ueIdx)->bool
  {
    v1 = uEdges(ueIdx,0);
    v2 = uEdges(ueIdx,1);
    return true;
  };


  // lambda——post_collapse()， update quadric
  post_collapse = [&v1,&v2,&quadrics](
      const Eigen::MatrixXd &                             ,   /*vers*/
      const Eigen::MatrixXi &                             ,   /*tris*/
      const Eigen::MatrixXi &                             ,   /*uEdges*/
      const Eigen::VectorXi &                             ,/*edgeUeInfo*/
      const Eigen::MatrixXi &                             ,  /*EF*/
      const Eigen::MatrixXi &                             ,  /*EI*/
      const igl::min_heap< std::tuple<double,int,int> > & ,/*Q*/
      const Eigen::VectorXi &                             ,/*EQ*/
      const Eigen::MatrixXd &                             ,   /*C*/
      const int                                           ,   /*ueIdx*/
      const int                                           ,  /*e1*/
      const int                                           ,  /*e2*/
      const int                                           ,  /*f1*/
      const int                                           ,  /*f2*/
      const bool                                          collapsed
      )->void
  {
    if(collapsed)
          quadrics[v1<v2?v1:v2] = quadrics[v1] + quadrics[v2];
  };

}

