// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_AVERAGE_ONTO_FACES_H
#define IGL_AVERAGE_ONTO_FACES_H
#include "igl_inline.h"

#include <Eigen/Dense>
namespace igl 
{
  // average_onto_vertices 
  // Move a scalar field defined on faces to vertices by averaging
  //
  // Input:
  //   tris  #tris by ss list of simples/faces
  //   S  #vers by dim list of per-vertex values
  // Output:
  //   SF  #tris by dim list of per-face values
  template <typename DerivedF, typename DerivedS, typename DerivedSF>
  IGL_INLINE void average_onto_faces(
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedS> & S,
    Eigen::PlainObjectBase<DerivedSF> & SF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "average_onto_faces.cpp"
#endif

#endif
