// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "boundary_loop.h"
#include "slice.h"
#include "triangle_triangle_adjacency.h"
#include "vertex_triangle_adjacency.h"
#include "is_border_vertex.h"
#include <set>

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
    const Eigen::MatrixBase<DerivedF> & tris,
    std::vector<std::vector<Index> >& L)
{
  using namespace std;
  using namespace Eigen;

  if(tris.rows() == 0)
    return;

  VectorXd Vdummy(tris.maxCoeff()+1,1);
  Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, Eigen::Dynamic> TT,TTi;
  vector<std::vector<int> > VF, VFi;
  triangle_triangle_adjacency(tris,TT,TTi);
  vertex_triangle_adjacency(Vdummy,tris,VF,VFi);

  vector<bool> unvisited = is_border_vertex(tris);
  set<int> unseen;
  for (size_t i = 0; i < unvisited.size(); ++i)
  {
    if (unvisited[i])
      unseen.insert(unseen.end(),i);
  }

  while (!unseen.empty())
  {
    vector<Index> l;

    // Get first vertex of loop
    int start = *unseen.begin();
    unseen.erase(unseen.begin());
    unvisited[start] = false;
    l.push_back(start);

    bool done = false;
    while (!done)
    {
      // Find next vertex
      bool newBndEdge = false;
      int v = l[l.size()-1];
      int next;
      for (int i = 0; i < (int)VF[v].size() && !newBndEdge; i++)
      {
        int fid = VF[v][i];

        if (TT.row(fid).minCoeff() < 0.) // Face contains boundary edge
        {
          int vLoc = -1;
          if (tris(fid,0) == v) vLoc = 0;
          if (tris(fid,1) == v) vLoc = 1;
          if (tris(fid,2) == v) vLoc = 2;

          int vNext = tris(fid,(vLoc + 1) % tris.cols());

          newBndEdge = false;
          if (unvisited[vNext] && TT(fid,vLoc) < 0)
          {
            next = vNext;
            newBndEdge = true;
          }
        }
      }

      if (newBndEdge)
      {
        l.push_back(next);
        unseen.erase(next);
        unvisited[next] = false;
      }
      else
        done = true;
    }
    L.push_back(l);
  }
}

template <typename DerivedF, typename Index>
IGL_INLINE void igl::boundary_loop(
  const Eigen::MatrixBase<DerivedF>& tris,
  std::vector<Index>& L)
{
  using namespace Eigen;
  using namespace std;

  if(tris.rows() == 0)
    return;

  vector<vector<int> > Lall;
  boundary_loop(tris,Lall);

  int idxMax = -1;
  size_t maxLen = 0;
  for (size_t i = 0; i < Lall.size(); ++i)
  {
    if (Lall[i].size() > maxLen)
    {
      maxLen = Lall[i].size();
      idxMax = i;
    }
  }

  //Check for meshes without boundary
  if (idxMax == -1)
  {
      L.clear();
      return;
  }

  L.resize(Lall[idxMax].size());
  for (size_t i = 0; i < Lall[idxMax].size(); ++i)
  {
    L[i] = Lall[idxMax][i];
  }
}

template <typename DerivedF, typename DerivedL>
IGL_INLINE void igl::boundary_loop(
  const Eigen::MatrixBase<DerivedF>& tris,
  Eigen::PlainObjectBase<DerivedL>& L)
{
  using namespace Eigen;
  using namespace std;

  if(tris.rows() == 0)
    return;

  vector<int> Lvec;
  boundary_loop(tris,Lvec);

  L.resize(Lvec.size(), 1);
  for (size_t i = 0; i < Lvec.size(); ++i)
    L(i) = Lvec[i];
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::boundary_loop<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::boundary_loop<Eigen::Matrix<int, -1, -1, 0, -1, -1>, int>(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&);
#endif
