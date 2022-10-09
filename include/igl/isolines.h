// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef IGL_ISOLINES_H
#define IGL_ISOLINES_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace igl
{
    // Constructs isolines for a function z given on a mesh (vers,tris)
    //
    //
    // Inputs:
    //   vers  #vers by dim list of mesh vertex positions
    //   tris  #tris by 3 list of mesh faces (must be triangles)
    //   z  #vers by 1 list of function values evaluated at vertices
    //   n  the number of desired isolines
    // Outputs:
    //   isoV  #isoV by dim list of isoline vertex positions
    //   isoE  #isoE by 2 list of isoline edge positions
    //
    //
    
    template <typename DerivedV,
    typename DerivedF,
    typename DerivedZ,
    typename DerivedIsoV,
    typename DerivedIsoE>
    IGL_INLINE void isolines(
                             const Eigen::MatrixBase<DerivedV>& vers,
                             const Eigen::MatrixBase<DerivedF>& tris,
                             const Eigen::MatrixBase<DerivedZ>& z,
                             const int n,
                             Eigen::PlainObjectBase<DerivedIsoV>& isoV,
                             Eigen::PlainObjectBase<DerivedIsoE>& isoE);
}

#ifndef IGL_STATIC_LIBRARY
#  include "isolines.cpp"
#endif

#endif
