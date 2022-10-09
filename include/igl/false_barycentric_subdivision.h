// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ADD_BARYCENTER_H
#define IGL_ADD_BARYCENTER_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl
{
  // Refine the mesh by adding the barycenter of each face
  // Inputs:
  //   vers       #vers by 3 coordinates of the vertices
  //   tris       #tris by 3 list of mesh faces (must be triangles)
  // Outputs:
  //   VD      #vers + #tris by 3 coordinate of the vertices of the dual mesh
  //           The added vertices are added at the end of VD (should not be
  //           same references as (vers,tris)
  //   FD      #tris*3 by 3 faces of the dual mesh
  //
  template <typename Scalar, typename Index>
  IGL_INLINE void false_barycentric_subdivision(
    const Eigen::PlainObjectBase<Scalar> & vers,
    const Eigen::PlainObjectBase<Index> & tris,
    Eigen::PlainObjectBase<Scalar> & VD,
    Eigen::PlainObjectBase<Index> & FD);

}

#ifndef IGL_STATIC_LIBRARY
#  include "false_barycentric_subdivision.cpp"
#endif

#endif
