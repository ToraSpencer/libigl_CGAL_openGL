// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LBS_MATRIX_H
#define IGL_LBS_MATRIX_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // LBS_MATRIX Linear blend skinning can be expressed by vers' = M * T where vers' is
  // a #vers by dim matrix of deformed vertex positions (one vertex per row), M is a
  // #vers by (dim+1)*#T (composed of weights and rest positions) and T is a
  // #T*(dim+1) by dim matrix of #T stacked transposed transformation matrices.
  // See equations (1) and (2) in "Fast Automatic Skinning Transformations"
  // [Jacobson et al 2012]
  //
  // Inputs:
  //   vers  #vers by dim list of rest positions
  //   W  #vers+ by #T  list of weights
  // Outputs:
  //   M  #vers by #T*(dim+1)
  //
  // In MATLAB:
  //   kron(ones(1,size(W,2)),[vers ones(size(vers,1),1)]).*kron(W,ones(1,size(vers,2)+1))
  IGL_INLINE void lbs_matrix(
    const Eigen::MatrixXd & vers, 
    const Eigen::MatrixXd & W,
    Eigen::MatrixXd & M);
  // LBS_MATRIX  construct a matrix that when multiplied against a column of
  // affine transformation entries computes new coordinates of the vertices
  //
  // I'm not sure it makes since that the result is stored as a sparse matrix.
  // The number of non-zeros per row *is* dependent on the number of mesh
  // vertices and handles.
  //
  // Inputs:
  //   vers  #vers by dim list of vertex rest positions
  //   W  #vers by #handles list of correspondence weights
  // Output:
  //   M  #vers * dim by #handles * dim * (dim+1) matrix such that
  //     new_V(:) = LBS(vers,W,A) = reshape(M * A,size(vers)), where A is a column
  //     vectors formed by the entries in each handle's dim by dim+1 
  //     transformation matrix. Specifcally, A =
  //       reshape(permute(Astack,[3 1 2]),n*dim*(dim+1),1)
  //     or A = [Lxx;Lyx;Lxy;Lyy;tx;ty], and likewise for other dim
  //     if Astack(:,:,i) is the dim by (dim+1) transformation at handle i
  IGL_INLINE void lbs_matrix_column(
    const Eigen::MatrixXd & vers, 
    const Eigen::MatrixXd & W,
    Eigen::SparseMatrix<double>& M);
  IGL_INLINE void lbs_matrix_column(
    const Eigen::MatrixXd & vers, 
    const Eigen::MatrixXd & W,
    Eigen::MatrixXd & M);
  // Same as LBS_MATRIX above but instead of giving W as a full matrix of weights
  // (each vertex has #handles weights), a constant number of weights are given
  // for each vertex.
  // 
  // Inputs:
  //   vers  #vers by dim list of vertex rest positions
  //   W  #vers by k  list of k correspondence weights per vertex
  //   WI  #vers by k  list of k correspondence weight indices per vertex. Such that
  //     W(j,WI(i)) gives the ith most significant correspondence weight on vertex j
  // Output:
  //   M  #vers * dim by #handles * dim * (dim+1) matrix such that
  //     new_V(:) = LBS(vers,W,A) = reshape(M * A,size(vers)), where A is a column
  //     vectors formed by the entries in each handle's dim by dim+1 
  //     transformation matrix. Specifcally, A =
  //       reshape(permute(Astack,[3 1 2]),n*dim*(dim+1),1)
  //     or A = [Lxx;Lyx;Lxy;Lyy;tx;ty], and likewise for other dim
  //     if Astack(:,:,i) is the dim by (dim+1) transformation at handle i
  //
  IGL_INLINE void lbs_matrix_column(
    const Eigen::MatrixXd & vers, 
    const Eigen::MatrixXd & W,
    const Eigen::MatrixXi & WI,
    Eigen::SparseMatrix<double>& M);
  IGL_INLINE void lbs_matrix_column(
    const Eigen::MatrixXd & vers, 
    const Eigen::MatrixXd & W,
    const Eigen::MatrixXi & WI,
    Eigen::MatrixXd & M);
}
#ifndef IGL_STATIC_LIBRARY
#include "lbs_matrix.cpp"
#endif
#endif
