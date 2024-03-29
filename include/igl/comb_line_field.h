// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Nico Pietroni <nico.pietroni@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_COMB_LINE_FIELD_H
#define IGL_COMB_LINE_FIELD_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Computes principal matchings of the vectors of a cross field across face edges,
  // and generates a combed cross field defined on the mesh faces

  // Inputs:
  //   vers          #vers by 3 eigen Matrix of mesh vertex 3D positions
  //   tris          #tris by 4 eigen Matrix of face (quad) indices
  //   PD1in      #tris by 3 eigen Matrix of the first per face cross field vector
  // Output:
  //   PD1out      #tris by 3 eigen Matrix of the first combed cross field vector


  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void comb_line_field(const Eigen::MatrixBase<DerivedV> &vers,
                                  const Eigen::MatrixBase<DerivedF> &tris,
                                  const Eigen::MatrixBase<DerivedV> &PD1in,
                                  Eigen::PlainObjectBase<DerivedV> &PD1out);
}
#ifndef IGL_STATIC_LIBRARY
#include "comb_line_field.cpp"
#endif

#endif
