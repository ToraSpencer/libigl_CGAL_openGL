// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MASSMATRIX_INTRINSIC_H
#define IGL_MASSMATRIX_INTRINSIC_H
#include "igl_inline.h"
#include "massmatrix.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{

  // Constructs the mass (area) matrix for a given mesh (vers,tris).
  //
  // Inputs:
  //   l  #l by simplex_size list of mesh edge lengths
  //   tris  #tris by simplex_size list of mesh elements (triangles or tetrahedra)
  //   type  one of the following ints:
  //     MASSMATRIX_TYPE_BARYCENTRIC  barycentric
  //     MASSMATRIX_TYPE_VORONOI voronoi-hybrid {default}
  //     MASSMATRIX_TYPE_FULL full {not implemented}
  // Outputs: 
  //   M  #vers by #vers mass matrix
  //
  // See also: adjacency_matrix
  //
  template <typename Derivedl, typename DerivedF, typename Scalar>
  IGL_INLINE void massmatrix_intrinsic(
    const Eigen::MatrixBase<Derivedl> & l, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    const MassMatrixType type,
    Eigen::SparseMatrix<Scalar>& M);
  // Inputs:
  //   n  number of vertices (>= tris.maxCoeff()+1)
  template <typename Derivedl, typename DerivedF, typename Scalar>
  IGL_INLINE void massmatrix_intrinsic(
    const Eigen::MatrixBase<Derivedl> & l, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    const MassMatrixType type,
    const int n,
    Eigen::SparseMatrix<Scalar>& M);
}

#ifndef IGL_STATIC_LIBRARY
#  include "massmatrix_intrinsic.cpp"
#endif

#endif


