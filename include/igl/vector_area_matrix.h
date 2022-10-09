// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_VECTOR_AREA_MATRIX_H
#define IGL_VECTOR_AREA_MATRIX_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // Constructs the symmetric area matrix A, s.t.  [vers.col(0)' vers.col(1)'] * A *
  // [vers.col(0); vers.col(1)] is the **vector area** of the mesh (vers,tris).
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for vers (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for tris (e.g. derived from
  //     MatrixXi)
  //   Scalar  scalar type for eigen sparse matrix (e.g. double)
  // Inputs:
  //   tris  #tris by 3 list of mesh faces (must be triangles)
  // Outputs:
  //   A  #Vx2 by #Vx2 area matrix
  //
  template <typename DerivedF, typename Scalar>
  IGL_INLINE void vector_area_matrix(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::SparseMatrix<Scalar>& A);
}

#ifndef IGL_STATIC_LIBRARY
#  include "vector_area_matrix.cpp"
#endif

#endif
