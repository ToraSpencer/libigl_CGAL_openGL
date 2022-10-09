// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ORIENT_OUTWARD_H
#define IGL_ORIENT_OUTWARD_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Orient each component (identified by C) of a mesh (vers,tris) so the normals on
  // average point away from the patch's centroid.
  //
  // Inputs:
  //   vers  #vers by 3 list of vertex positions
  //   tris  #tris by 3 list of triangle indices
  //   C  #tris list of components (output of orientable_patches)
  // Outputs:
  //   FF  #tris by 3 list of new triangle indices such that FF(~I,:) = tris(~I,:) and
  //     FF(I,:) = fliplr(tris(I,:)) (OK if &FF = &tris)
  //   I  max(C)+1 list of whether face has been flipped
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedC,
    typename DerivedFF,
    typename DerivedI>
  IGL_INLINE void orient_outward(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedI> & I);
};

#ifndef IGL_STATIC_LIBRARY
#  include "orient_outward.cpp"
#endif

#endif
