// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ALL_PAIRS_DISTANCES_H
#define IGL_ALL_PAIRS_DISTANCES_H
#include "igl_inline.h"

namespace igl
{
  // ALL_PAIRS_DISTANCES compute distances between each point i in vers and point j
  // in U
  // 
  // D = all_pairs_distances(vers,U)
  // 
  // Templates:
  //   Mat  matrix class like MatrixXd
  // Inputs:
  //   vers  #vers by dim list of points
  //   U  #U by dim list of points
  //   squared  whether to return squared distances
  // Outputs:
  //   D  #vers by #U matrix of distances, where D(i,j) gives the distance or
  //     squareed distance between vers(i,:) and U(j,:)
  // 
  template <typename Mat>
  IGL_INLINE void all_pairs_distances(
    const Mat & vers,
    const Mat & U,
    const bool squared, 
    Mat & D);
}

#ifndef IGL_STATIC_LIBRARY
#  include "all_pairs_distances.cpp"
#endif

#endif
