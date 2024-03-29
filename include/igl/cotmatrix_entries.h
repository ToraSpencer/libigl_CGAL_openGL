// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COTMATRIX_ENTRIES_H
#define IGL_COTMATRIX_ENTRIES_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // COTMATRIX_ENTRIES compute the cotangents of each angle in mesh (vers,tris)
  // 
  // Inputs:
  //   vers  #vers by dim list of rest domain positions
  //   tris  #tris by {3|4} list of {triangle|tetrahedra} indices into vers
  // Outputs:
  //     C  #tris by 3 list of 1/2*cotangents corresponding angles
  //       for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //   OR
  //     C  #tris by 6 list of 1/6*cotangents of dihedral angles*edge lengths
  //       for tets, columns along edges [1,2],[2,0],[0,1],[3,0],[3,1],[3,2] 
  //
  template <typename DerivedV, typename DerivedF, typename DerivedC>
  IGL_INLINE void cotmatrix_entries(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedC>& C);
  // Intrinsic version.
  //
  // Inputs:
  //   l  #tris by 3 list of triangle edge lengths (see edge_lengths)
  // Outputs:
  //   C  #tris by 3 list of 1/2*cotangents corresponding angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  template <typename Derivedl, typename DerivedC>
  IGL_INLINE void cotmatrix_entries(
    const Eigen::MatrixBase<Derivedl>& l,
    Eigen::PlainObjectBase<DerivedC>& C);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cotmatrix_entries.cpp"
#endif

#endif
