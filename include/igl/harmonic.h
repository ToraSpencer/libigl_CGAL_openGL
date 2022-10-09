// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_HARMONIC_H
#define IGL_HARMONIC_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // Compute k-harmonic weight functions "coordinates".
  //
  //
  // Inputs:
  //   vers  #vers by dim vertex positions
  //   tris  #tris by simplex-size list of element indices
  //   b  #b boundary indices into vers
  //   bc #b by #W list of boundary values
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   W  #vers by #W list of weights
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedb,
    typename Derivedbc,
    typename DerivedW>
  IGL_INLINE bool harmonic(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<Derivedb> & b,
    const Eigen::MatrixBase<Derivedbc> & bc,
    const int k,
    Eigen::PlainObjectBase<DerivedW> & W);
  // Compute harmonic map using uniform laplacian operator
  //
  // Inputs:
  //   tris  #tris by simplex-size list of element indices
  //   b  #b boundary indices into vers
  //   bc #b by #W list of boundary values
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   W  #vers by #W list of weights
  //
  template <
    typename DerivedF,
    typename Derivedb,
    typename Derivedbc,
    typename DerivedW>
  IGL_INLINE bool harmonic(
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<Derivedb> & b,
    const Eigen::MatrixBase<Derivedbc> & bc,
    const int k,
    Eigen::PlainObjectBase<DerivedW> & W);
  // Compute a harmonic map using a given Laplacian and mass matrix
  //
  // Inputs:
  //   L  #vers by #vers discrete (integrated) Laplacian  
  //   M  #vers by #vers mass matrix
  //   b  #b boundary indices into vers
  //   bc  #b by #W list of boundary values
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   W  #vers by #vers list of weights
  template <
    typename DerivedL,
    typename DerivedM,
    typename Derivedb,
    typename Derivedbc,
    typename DerivedW>
  IGL_INLINE bool harmonic(
    const Eigen::SparseCompressedBase<DerivedL> & L,
    const Eigen::SparseCompressedBase<DerivedM> & M,
    const Eigen::MatrixBase<Derivedb> & b,
    const Eigen::MatrixBase<Derivedbc> & bc,
    const int k,
    Eigen::PlainObjectBase<DerivedW> & W);
  // Build the discrete k-harmonic operator (computing integrated quantities).
  // That is, if the k-harmonic PDE is Q x = 0, then this minimizes x' Q x
  //
  // Inputs:
  //   L  #vers by #vers discrete (integrated) Laplacian  
  //   M  #vers by #vers mass matrix
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   Q  #vers by #vers discrete (integrated) k-Laplacian  
  template <
    typename DerivedL,
    typename DerivedM,
    typename DerivedQ>
  IGL_INLINE void harmonic(
    const Eigen::SparseCompressedBase<DerivedL> & L,
    const Eigen::SparseCompressedBase<DerivedM> & M,
    const int k,
    DerivedQ & Q);
  // Inputs:
  //   vers  #vers by dim vertex positions
  //   tris  #tris by simplex-size list of element indices
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   Q  #vers by #vers discrete (integrated) k-Laplacian  
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedQ>
  IGL_INLINE void harmonic(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const int k,
    DerivedQ & Q);
};

#ifndef IGL_STATIC_LIBRARY
#include "harmonic.cpp"
#endif
#endif
