// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ARAP_RHS_H
#define IGL_ARAP_RHS_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/ARAPEnergyType.h>

namespace igl
{
  // ARAP_RHS build right-hand side constructor of global poisson solve for
  // various Arap energies
  // Inputs:
  //   vers  #vers by Vdim list of initial domain positions
  //   tris  #tris by 3 list of triangle indices into vers
  //   dim  dimension being used at solve time. For deformation usually dim =
  //     vers.cols(), for surface parameterization vers.cols() = 3 and dim = 2
  //   energy  igl::ARAPEnergyType enum value defining which energy is being
  //     used. See igl::ARAPEnergyType.h for valid options and explanations.
  // Outputs:
  //   K  #vers*dim by #(tris|vers)*dim*dim matrix such that:
  //     b = K * reshape(permute(R,[3 1 2]),size(vers|tris,1)*size(vers,2)*size(vers,2),1);
  //
  // See also: arap_linear_block
  template<typename DerivedV, typename DerivedF, typename DerivedK>
  IGL_INLINE void arap_rhs(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const int dim,
    const igl::ARAPEnergyType energy,
    Eigen::SparseCompressedBase<DerivedK>& K);
}
#ifndef IGL_STATIC_LIBRARY
#include "arap_rhs.cpp"
#endif
#endif
