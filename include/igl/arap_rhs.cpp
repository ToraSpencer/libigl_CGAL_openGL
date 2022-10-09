// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "arap_rhs.h"
#include "arap_linear_block.h"
#include "verbose.h"
#include "repdiag.h"
#include "cat.h"
#include <iostream>

template<typename DerivedV, typename DerivedF, typename DerivedK>
IGL_INLINE void igl::arap_rhs(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const int dim,
    const igl::ARAPEnergyType energy,
    Eigen::SparseCompressedBase<DerivedK>& K)
{
  using namespace std;
  using namespace Eigen;
  // Number of dimensions
  int Vdim = vers.cols();
  //// Number of mesh vertices
  //int n = vers.rows();
  //// Number of mesh elements
  //int m = tris.rows();
  //// number of rotations
  //int nr;
  switch(energy)
  {
    case ARAP_ENERGY_TYPE_SPOKES:
      //nr = n;
      break;
    case ARAP_ENERGY_TYPE_SPOKES_AND_RIMS:
      //nr = n;
      break;
    case ARAP_ENERGY_TYPE_ELEMENTS:
      //nr = m;
      break;
    default:
      fprintf(
        stderr,
        "arap_rhs.h: Error: Unsupported arap energy %d\n",
        energy);
      return;
  }

  DerivedK KX,KY,KZ;
  arap_linear_block(vers,tris,0,energy,KX);
  arap_linear_block(vers,tris,1,energy,KY);
  if(Vdim == 2)
  {
    K = cat(2,repdiag(KX,dim),repdiag(KY,dim));
  }else if(Vdim == 3)
  {
    arap_linear_block(vers,tris,2,energy,KZ);
    if(dim == 3)
    {
      K = cat(2,cat(2,repdiag(KX,dim),repdiag(KY,dim)),repdiag(KZ,dim));
    }else if(dim ==2)
    {
      DerivedK ZZ(KX.rows()*2,KX.cols());
      K = cat(2,cat(2,
            cat(2,repdiag(KX,dim),ZZ),
            cat(2,repdiag(KY,dim),ZZ)),
            cat(2,repdiag(KZ,dim),ZZ));
    }else
    {
      assert(false);
      fprintf(
      stderr,
      "arap_rhs.h: Error: Unsupported dimension %d\n",
      dim);
    }
  }else
  {
    assert(false);
    fprintf(
     stderr,
     "arap_rhs.h: Error: Unsupported dimension %d\n",
     Vdim);
    return;
  }

}



#ifdef IGL_STATIC_LIBRARY
template void igl::arap_rhs(const Eigen::MatrixBase<Eigen::MatrixXd> & vers, const Eigen::MatrixBase<Eigen::MatrixXi> & tris,const int dim, const igl::ARAPEnergyType energy,Eigen::SparseCompressedBase<Eigen::SparseMatrix<double>>& K);
#endif