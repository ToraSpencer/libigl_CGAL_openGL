// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "average_onto_vertices.h"

template<typename DerivedV,typename DerivedF,typename DerivedS,typename DerivedSV >
IGL_INLINE void igl::average_onto_vertices(const Eigen::MatrixBase<DerivedV> &vers,
  const Eigen::MatrixBase<DerivedF> &tris,
  const Eigen::MatrixBase<DerivedS> &S,
  Eigen::PlainObjectBase<DerivedSV> &SV)
{
  SV = DerivedS::Zero(vers.rows(),S.cols());
  Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> COUNT(vers.rows());
  COUNT.setZero();
  for (int i = 0; i <tris.rows(); ++i)
  {
    for (int j = 0; j<tris.cols(); ++j)
    {
      SV.row(tris(i,j)) += S.row(i);
      COUNT[tris(i,j)] ++;
    }
  }
  for (int i = 0; i <vers.rows(); ++i)
    SV.row(i) /= COUNT[i];
};

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
