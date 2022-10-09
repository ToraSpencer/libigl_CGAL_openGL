// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "false_barycentric_subdivision.h"

#include "verbose.h"
#include <algorithm>
#include <igl/barycenter.h>

template <typename Scalar, typename Index>
IGL_INLINE void igl::false_barycentric_subdivision(
    const Eigen::PlainObjectBase<Scalar> & vers,
    const Eigen::PlainObjectBase<Index> & tris,
    Eigen::PlainObjectBase<Scalar> & VD,
    Eigen::PlainObjectBase<Index> & FD)
{
  using namespace Eigen;
  // Compute face barycenter
  Eigen::MatrixXd BC;
  igl::barycenter(vers,tris,BC);

  // Add the barycenters to the vertices
  VD.resize(vers.rows()+tris.rows(),3);
  VD.block(0,0,vers.rows(),3) = vers;
  VD.block(vers.rows(),0,tris.rows(),3) = BC;

  // Each face is split four ways
  FD.resize(tris.rows()*3,3);

  for (unsigned i=0; i<tris.rows(); ++i)
  {
    int i0 = tris(i,0);
    int i1 = tris(i,1);
    int i2 = tris(i,2);
    int i3 = vers.rows() + i;

    Vector3i F0,F1,F2;
    F0 << i0,i1,i3;
    F1 << i1,i2,i3;
    F2 << i2,i0,i3;

    FD.row(i*3 + 0) = F0;
    FD.row(i*3 + 1) = F1;
    FD.row(i*3 + 2) = F2;
  }


}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::false_barycentric_subdivision<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template void igl::false_barycentric_subdivision<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
