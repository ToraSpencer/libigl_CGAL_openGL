// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "faces_first.h"

#include <vector>
#include <Eigen/Dense>

template <typename MatV, typename MatF, typename VecI>
IGL_INLINE void igl::faces_first(
  const MatV & vers, 
  const MatF & tris, 
  MatV & RV, 
  MatF & RF, 
  VecI & IM)
{
  assert(&vers != &RV);
  assert(&tris != &RF);
  using namespace std;
  using namespace Eigen;
  vector<bool> in_face(vers.rows());
  for(int i = 0; i<tris.rows(); i++)
  {
    for(int j = 0; j<tris.cols(); j++)
    {
      in_face[tris(i,j)] = true;
    }
  }
  // count number of vertices not in faces
  int num_in_F = 0;
  for(int i = 0;i<vers.rows();i++)
  {
    num_in_F += (in_face[i]?1:0);
  }
  // list of unique vertices that occur in tris
  VectorXi U(num_in_F);
  // list of unique vertices that do not occur in tris
  VectorXi NU(vers.rows()-num_in_F);
  int Ui = 0;
  int NUi = 0;
  // loop over vertices
  for(int i = 0;i<vers.rows();i++)
  {
    if(in_face[i])
    {
      U(Ui) = i;
      Ui++;
    }else
    {
      NU(NUi) = i;
      NUi++;
    }
  }
  IM.resize(vers.rows());
  // reindex vertices that occur in faces to be first
  for(int i = 0;i<U.size();i++)
  {
    IM(U(i)) = i;
  }
  // reindex vertices that do not occur in faces to come after those that do
  for(int i = 0;i<NU.size();i++)
  {
    IM(NU(i)) = i+U.size();
  }
  RF.resizeLike(tris);
  // Reindex faces
  for(int i = 0; i<tris.rows(); i++)
  {
    for(int j = 0; j<tris.cols(); j++)
    {
      RF(i,j) = IM(tris(i,j));
    }
  }
  RV.resizeLike(vers);
  // Reorder vertices
  for(int i = 0;i<vers.rows();i++)
  {
    RV.row(IM(i)) = vers.row(i);
  }
}

template <typename MatV, typename MatF, typename VecI>
IGL_INLINE void igl::faces_first(
  MatV & vers, 
  MatF & tris, 
  VecI & IM)
{
  MatV RV;
  // Copying tris may not be needed, seems RF = tris is safe (whereas RV = vers is not)
  MatF RF;
  igl::faces_first(vers,tris,RV,RF,IM);
  vers = RV;
  tris = RF;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::faces_first<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::Matrix<double, -1, -1, 0, -1, -1>&, Eigen::Matrix<int, -1, -1, 0, -1, -1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);
#endif
