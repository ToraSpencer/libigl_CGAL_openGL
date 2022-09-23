#ifndef IGL_SIMPLIFY_POLYHEDRON_H
#define IGL_SIMPLIFY_POLYHEDRON_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Simplify a polyhedron represented as a triangle mesh (OV,OF) by collapsing
  // any edge that doesn't contribute to defining surface's pointset. This
  // _would_ also make sense for open and non-manifold meshes, but the current
  // implementation only works with closed manifold surfaces with well defined
  // triangle normals.
  //
  // Inputs:
  //   OV  #OV by 3 list of input mesh vertex positions
  //   OF  #OF by 3 list of input mesh triangle indices into OV
  // Outputs:
  //   V  #V by 3 list of output mesh vertex positions
  //   F  #F by 3 list of input mesh triangle indices into V
  //   J  #F list of indices into OF of birth parents
  IGL_INLINE void simplify_polyhedron(
    const Eigen::MatrixXd & OV,
    const Eigen::MatrixXi & OF,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::VectorXi & J);
}
#ifndef IGL_STATIC_LIBRARY
#  include "simplify_polyhedron.cpp"
#endif
#endif
