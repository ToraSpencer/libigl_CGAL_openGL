// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "orient_outward.h"
#include "per_face_normals.h"
#include "barycenter.h"
#include "doublearea.h"
#include <iostream>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedFF,
  typename DerivedI>
IGL_INLINE void igl::orient_outward(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  const Eigen::MatrixBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace std;
  assert(C.rows() == tris.rows());
  assert(tris.cols() == 3);
  assert(vers.cols() == 3);

  // number of faces
  const int m = tris.rows();
  // number of patches
  const int num_cc = C.maxCoeff()+1;
  I.resize(num_cc);
  if(&FF != &tris)
  {
    FF = tris;
  }
  DerivedV N,BC,BCmean;
  Matrix<typename DerivedV::Scalar,Dynamic,1> A;
  VectorXd totA(num_cc), dot(num_cc);
  Matrix<typename DerivedV::Scalar,3,1> Z(1,1,1);
  per_face_normals(vers,tris,Z.normalized(),N);
  barycenter(vers,tris,BC);
  doublearea(vers,tris,A);
  BCmean.setConstant(num_cc,3,0);
  dot.setConstant(num_cc,1,0);
  totA.setConstant(num_cc,1,0);
  // loop over faces
  for(int f = 0;f<m;f++)
  {
    BCmean.row(C(f)) += A(f)*BC.row(f);
    totA(C(f))+=A(f);
  }
  // take area weighted average
  for(int c = 0;c<num_cc;c++)
  {
    BCmean.row(c) /= (typename DerivedV::Scalar) totA(c);
  }
  // subtract bcmean
  for(int f = 0;f<m;f++)
  {
    BC.row(f) -= BCmean.row(C(f));
    dot(C(f)) += A(f)*N.row(f).dot(BC.row(f));
  }
  // take area weighted average
  for(int c = 0;c<num_cc;c++)
  {
    dot(c) /= (typename DerivedV::Scalar) totA(c);
    if(dot(c) < 0)
    {
      I(c) = true;
    }else
    {
      I(c) = false;
    }
  }
  // flip according to I
  for(int f = 0;f<m;f++)
  {
    if(I(C(f)))
    {
      FF.row(f) = FF.row(f).reverse().eval();
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::orient_outward<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::orient_outward<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif

