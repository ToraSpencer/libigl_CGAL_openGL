// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sample_edges.h"

IGL_INLINE void igl::sample_edges(
  const Eigen::MatrixXd & vers,
  const Eigen::MatrixXi & E,
  const int k,
  Eigen::MatrixXd & S)
{
  using namespace Eigen;
  // Resize output
  S.resize(vers.rows() + k * E.rows(),vers.cols());
  // Copy vers at front of S
  S.block(0,0,vers.rows(),vers.cols()) = vers;

  // loop over edges
  for(int i = 0;i<E.rows();i++)
  {
    VectorXd tip = vers.row(E(i,0));
    VectorXd tail = vers.row(E(i,1));
    for(int s=0;s<k;s++)
    {
      double f = double(s+1)/double(k+1);
      S.row(vers.rows()+k*i+s) = f*tail + (1.0-f)*tip;
    }
  }
}
