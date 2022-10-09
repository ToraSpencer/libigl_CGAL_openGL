// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FACET_ADJACENCY_MATRIX_H
#define IGL_FACET_ADJACENCY_MATRIX_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "igl_inline.h"

namespace igl
{
  // Construct a #tris×#tris adjacency matrix with A(i,j)>0 indicating that faces i and j
  // share an edge.
  //
  // Inputs:
  //   tris  #tris by 3 list of facets
  // Outputs:
  //   A  #tris by #tris adjacency matrix
  template <typename DerivedF, typename Atype>
  IGL_INLINE void facet_adjacency_matrix(
    const Eigen::MatrixBase<DerivedF> & tris,
    Eigen::SparseMatrix<Atype> & A);
};

#ifndef IGL_STATIC_LIBRARY
#  include "facet_adjacency_matrix.cpp"
#endif

#endif
