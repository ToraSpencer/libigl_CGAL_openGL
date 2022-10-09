// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "polar_dec.h"
#include "polar_svd.h"
#ifdef _WIN32
#else
#  include <fenv.h>
#endif
#include <cmath>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cfenv>

// From Olga's CGAL mentee's ARAP code
template <
  typename DerivedA,
  typename DerivedR,
  typename DerivedT,
  typename DerivedU,
  typename DerivedS,
  typename DerivedV>
IGL_INLINE void igl::polar_dec(
  const Eigen::PlainObjectBase<DerivedA> & A,
  Eigen::PlainObjectBase<DerivedR> & R,
  Eigen::PlainObjectBase<DerivedT> & T,
  Eigen::PlainObjectBase<DerivedU> & U,
  Eigen::PlainObjectBase<DerivedS> & S,
  Eigen::PlainObjectBase<DerivedV> & vers)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedA::Scalar Scalar;

  const Scalar th = std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());

  Eigen::SelfAdjointEigenSolver<DerivedA> eig;
  feclearexcept(FE_UNDERFLOW);
  eig.computeDirect(A.transpose()*A);
  if(fetestexcept(FE_UNDERFLOW) || eig.eigenvalues()(0)/eig.eigenvalues()(2)<th)
  {
    cout<<"resorting to svd 1..."<<endl;
    return polar_svd(A,R,T,U,S,vers);
  }

  S = eig.eigenvalues().cwiseSqrt();

  vers = eig.eigenvectors();
  U = A * vers;
  R = U * S.asDiagonal().inverse() * vers.transpose();
  T = vers * S.asDiagonal() * vers.transpose();

  S = S.reverse().eval();
  vers = vers.rowwise().reverse().eval();
  U = U.rowwise().reverse().eval() * S.asDiagonal().inverse();

  if(R.determinant() < 0)
  {
    // Annoyingly the .eval() is necessary
    auto W = vers.eval();
    const auto & SVT = S.asDiagonal() * vers.adjoint();
    W.col(vers.cols()-1) *= -1.;
    R = U*W.transpose();
    T = W*SVT;
  }

  if(std::fabs(R.squaredNorm()-3.) > th)
  {
    cout<<"resorting to svd 2..."<<endl;
    return polar_svd(A,R,T,U,S,vers);
  }
}

template <
  typename DerivedA,
  typename DerivedR,
  typename DerivedT>
IGL_INLINE void igl::polar_dec(
  const Eigen::PlainObjectBase<DerivedA> & A,
  Eigen::PlainObjectBase<DerivedR> & R,
  Eigen::PlainObjectBase<DerivedT> & T)
{
  DerivedA U;
  DerivedA vers;
  Eigen::Matrix<typename DerivedA::Scalar,DerivedA::RowsAtCompileTime,1> S;
  return igl::polar_dec(A,R,T,U,S,vers);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template  void igl::polar_dec<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::polar_dec<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 2, 2, 0, 2, 2>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 2, 0, 2, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
