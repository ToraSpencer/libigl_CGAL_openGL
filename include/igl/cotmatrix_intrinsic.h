// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COTMATRIX_INTRINSIC_H
#define IGL_COTMATRIX_INTRINSIC_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  // Constructs the cotangent stiffness matrix (discrete laplacian) for a given
  // mesh with faces tris and edge lengths l.
  //
  // Inputs:
  //   l  #tris by 3 list of (half-)edge lengths
  //   tris  #tris by 3 list of face indices into some (not necessarily
  //     determined/embedable) list of vertex positions vers. It is assumed #vers ==
  //     tris.maxCoeff()+1
  // Outputs:
  //   L  #vers by #vers sparse Laplacian matrix
  //
  // See also: cotmatrix, intrinsic_delaunay_cotmatrix
  template <typename Derivedl, typename DerivedF, typename Scalar>
  IGL_INLINE void cotmatrix_intrinsic(
    const Eigen::MatrixBase<Derivedl> & l, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    Eigen::SparseMatrix<Scalar>& L);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cotmatrix_intrinsic.cpp"
#endif

#endif
