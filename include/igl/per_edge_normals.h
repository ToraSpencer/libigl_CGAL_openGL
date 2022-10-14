// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PER_EDGE_NORMALS_H
#define IGL_PER_EDGE_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  enum PerEdgeNormalsWeightingType
  {
    // Incident face normals have uniform influence on edge normal
    PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM = 0,
    // Incident face normals are averaged weighted by area
    PER_EDGE_NORMALS_WEIGHTING_TYPE_AREA = 1,
    // Area weights
    PER_EDGE_NORMALS_WEIGHTING_TYPE_DEFAULT = 2,
    NUM_PER_EDGE_NORMALS_WEIGHTING_TYPE = 3
  };
  // Compute face normals via vertex position list, face list
  // Inputs:
  //   vers  #vers by 3 eigen Matrix of mesh vertex 3D positions
  //   tris  #tris by 3 eigen Matrix of face (triangle) indices
  //   weight  weighting type
  //   FN  #tris by 3 matrix of 3D face normals per face
  // Output:
  //   N  #2 by 3 matrix of mesh edge 3D normals per row
  //   E  #E by 2 matrix of edge indices per row
  //   edgeUeInfo  #E by 1 matrix of indices from all edges to E
  //
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedFN,
    typename DerivedN,
    typename DerivedE,
    typename DerivedEMAP>
  IGL_INLINE void per_edge_normals(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const PerEdgeNormalsWeightingType weight,
    const Eigen::MatrixBase<DerivedFN>& FN,
    Eigen::PlainObjectBase<DerivedN> & N,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo);
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedN,
    typename DerivedE,
    typename DerivedEMAP>
  IGL_INLINE void per_edge_normals(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const PerEdgeNormalsWeightingType weight,
    Eigen::PlainObjectBase<DerivedN> & N,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo);
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedN,
    typename DerivedE,
    typename DerivedEMAP>
  IGL_INLINE void per_edge_normals(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedN> & N,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo);
}

#ifndef IGL_STATIC_LIBRARY
#  include "per_edge_normals.cpp"
#endif

#endif
