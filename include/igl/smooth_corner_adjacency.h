#ifndef IGL_SMOOTH_CORNER_ADJACENCY_H
#define IGL_SMOOTH_CORNER_ADJACENCY_H

#include <Eigen/Core>
#include "igl_inline.h"

namespace igl
{
  // Determine the corner-to-face adjacency relationship that can be used for
  // computing crease-aware per-corner normals.
  //
  // Inputs:
  //   vers  #vers by 3 list of mesh vertex positions
  //   tris  #tris by 3 list of triangle mesh indices into rows of vers
  //   corner_threshold_radians  dihedral angle considered non-smooth (in
  //     radians)
  // Outputs:
  //   CI  #CI list of face neighbors as indices into rows of tris
  //   CC  3*#tris+1 list of cumulative sizes so that CC(i*3+j+1) - CC(i*3+j) is
  //     the number of faces considered smoothly incident on corner at tris(i,j)
  void smooth_corner_adjacency(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const double corner_threshold_radians,
    Eigen::VectorXi & CI,
    Eigen::VectorXi & CC);
  // Determine the effective corner-to-face adjacency relationship implied by a
  // set of indexed vertex positions (FV) and normals (FV) (e.g., those read in
  // from a .obj file).
  //
  // Inputs:
  //   FV  #tris by 3 list of triangle mesh indices into rows of some vers
  //   FN  #tris by 3 list of triangle mesh indices into rows of some N
  // Outputs:
  //   CI  #CI list of face neighbors as indices into rows of tris
  //   CC  3*#tris+1 list of cumulative sizes so that CC(i*3+j+1) - CC(i*3+j) is
  //     the number of faces considered smoothly incident on corner at tris(i,j)
  void smooth_corner_adjacency(
    const Eigen::MatrixXi & FV,
    const Eigen::MatrixXi & FN,
    Eigen::VectorXi & CI,
    Eigen::VectorXi & CC);
}

#ifndef IGL_STATIC_LIBRARY
#  include "smooth_corner_adjacency.cpp"
#endif

#endif
