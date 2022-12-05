#include "max_faces_stopping_condition.h"


IGL_INLINE void igl::max_faces_stopping_condition(
  int & trisCountNew,
  const int trisCount,
  const int tarTrisCount,
  decimate_stopping_condition_callback & stopping_condition)
{
  stopping_condition = 
    [trisCount,tarTrisCount,&trisCountNew](
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const igl::min_heap< std::tuple<double,int,int> > & ,
    const Eigen::VectorXi &                             ,
    const Eigen::MatrixXd &,
    const int,
    const int,
    const int,
    const int f1,
    const int f2)->bool
    {
      // Only subtract if we're collapsing a real face
      if(f1 < trisCount) 
          trisCountNew-=1;
      if(f2 < trisCount) 
          trisCountNew-=1;

      return trisCountNew<=(int)tarTrisCount;
    };
}


IGL_INLINE igl::decimate_stopping_condition_callback
igl::max_faces_stopping_condition(
  int & trisCountNew,
  const int trisCount,
  const int tarTrisCount)
{
  decimate_stopping_condition_callback stopping_condition;
  max_faces_stopping_condition(trisCountNew,trisCount,tarTrisCount,stopping_condition);
  return stopping_condition;
}
