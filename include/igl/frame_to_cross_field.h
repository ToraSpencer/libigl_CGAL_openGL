// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FRAME_TO_CROSS_FIELD_H
#define IGL_FRAME_TO_CROSS_FIELD_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl
{
  // Convert a frame field into its closest cross field
  // Inputs:
  //   vers       #vers by 3 coordinates of the vertices
  //   tris       #tris by 3 list of mesh faces (must be triangles)
  //   FF1     #tris by 3 the first representative vector of the frame field (up to permutation and sign)
  //   FF2     #tris by 3 the second representative vector of the frame field (up to permutation and sign)
  //
  // Outputs:
  //   X       #tris by 3 representative vector of the closest cross field
  //
  IGL_INLINE void frame_to_cross_field(
    const Eigen::MatrixXd& vers,
    const Eigen::MatrixXi& tris,
    const Eigen::MatrixXd& FF1,
    const Eigen::MatrixXd& FF2,
    Eigen::MatrixXd& X);

}

#ifndef IGL_STATIC_LIBRARY
#  include "frame_to_cross_field.cpp"
#endif

#endif
