// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITEMESH_H
#define IGL_WRITEMESH_H
#include "igl_inline.h"

#include <string>
#include <vector>
#include <Eigen/Core>

namespace igl
{
  // save a tetrahedral volume mesh to a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be cast as double)
  //   Index  type for indices (will be cast to int)
  // Input:
  //   mesh_file_name  path of .mesh file
  //   vers  double matrix of vertex positions  #vers by 3
  //   T  #T list of tet indices into vertex positions
  //   tris  #tris list of face indices into vertex positions
  //
  // Known bugs: Holes and regions are not supported
  template <typename Scalar, typename Index>
  IGL_INLINE bool writeMESH(
    const std::string mesh_file_name,
    const std::vector<std::vector<Scalar > > & vers,
    const std::vector<std::vector<Index > > & T,
    const std::vector<std::vector<Index > > & tris);

  // Templates:
  //   DerivedV  real-value: i.e. from MatrixXd
  //   DerivedT  integer-value: i.e. from MatrixXi
  //   DerivedF  integer-value: i.e. from MatrixXi
  // Input:
  //   mesh_file_name  path of .mesh file
  //   vers  eigen double matrix #vers by 3
  //   T  eigen int matrix #T by 4
  //   tris  eigen int matrix #tris by 3
  template <typename DerivedV, typename DerivedT, typename DerivedF>
  IGL_INLINE bool writeMESH(
    const std::string str,
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedT> & T,
    const Eigen::MatrixBase<DerivedF> & tris);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeMESH.cpp"
#endif

#endif
