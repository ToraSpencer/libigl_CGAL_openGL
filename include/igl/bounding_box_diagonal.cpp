// This file is part of libigl, a simple c++ geometry processing library.
#include "bounding_box_diagonal.h"
#include "mat_max.h"
#include "mat_min.h"
#include <cmath>

IGL_INLINE double igl::bounding_box_diagonal(
  const Eigen::MatrixXd & vers)
{
  using namespace Eigen;
  VectorXd maxV,minV;
  VectorXi maxVI,minVI;
  mat_max(vers,1,maxV,maxVI);
  mat_min(vers,1,minV,minVI);
  return sqrt((maxV-minV).array().square().sum());
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
