// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_PROGRESSIVE_HULLS_COST_AND_PLACEMENT_H
#define IGL_COPYLEFT_PROGRESSIVE_HULLS_COST_AND_PLACEMENT_H
#include <Eigen/Core>
#include "../igl_inline.h"
namespace igl
{
  namespace copyleft
  {
    // A "cost and placement" compatible with `igl::decimate` implementing the
    // "progressive hulls" algorithm in "Silhouette clipping" [Sander et al.
    // 2000]. This implementation fixes an issue that the original linear
    // program becomes unstable for flat patches by introducing a small
    // quadratic energy term pulling the collapsed edge toward its midpoint.
    // This function is not really meant to be called directly but rather
    // passed to `igl::decimate` as a handle.
    //
    // Inputs:
    //   e  index of edge to be collapsed
    //   vers  #vers by 3 list of vertex positions
    //   tris  #tris by 3 list of faces indices into vers
    //   E  #E by 3 list of edges indices into vers
    //   edgeUeInfo #tris*3 list of indices into E, mapping each directed edge to unique
    //     unique edge in E
    //   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
    //     tris(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
    //     e=(j->i)
    //   EI  #E by 2 list of edge flap corners (see above).
    // Outputs:
    //   cost  cost of collapsing edge e
    //   p  position to place collapsed vertex
    //
    IGL_INLINE void progressive_hulls_cost_and_placement(
      const int e,
      const Eigen::MatrixXd & vers,
      const Eigen::MatrixXi & tris,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & edgeUeInfo,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      double & cost,
      Eigen::RowVectorXd & p);
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "progressive_hulls_cost_and_placement.cpp"
#endif
#endif
