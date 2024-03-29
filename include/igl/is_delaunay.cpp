// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "is_delaunay.h"
#include "unique_edge_map.h"
#include <cassert>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedD>
IGL_INLINE void igl::is_delaunay(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  Eigen::PlainObjectBase<DerivedD> & D)
{
  typedef typename DerivedV::Scalar Scalar;
  // Should use Shewchuk's predicates instead.
  const auto float_incircle = [](
    const Scalar pa[2],
    const Scalar pb[2],
    const Scalar pc[2],
    const Scalar pd[2])->short
  {
    // I acknowledge that I am cating to double
    const Eigen::Matrix3d A = (Eigen::Matrix3d(3,3)<<
      pa[0]-pd[0], pa[1]-pd[1],(pa[0]-pd[0])*(pa[0]-pd[0])+(pa[1]-pd[1])*(pa[1]-pd[1]),
      pb[0]-pd[0], pb[1]-pd[1],(pb[0]-pd[0])*(pb[0]-pd[0])+(pb[1]-pd[1])*(pb[1]-pd[1]),
      pc[0]-pd[0], pc[1]-pd[1],(pc[0]-pd[0])*(pc[0]-pd[0])+(pc[1]-pd[1])*(pc[1]-pd[1])
      ).finished();
    const Scalar detA = A.determinant();
    return (Scalar(0) < detA) - (detA < Scalar(0));
  };

  typedef Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,2> MatrixX2I;
  typedef Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI edgeUeInfo;
  std::vector<std::vector<typename DerivedF::Scalar> > uE2E;
  igl::unique_edge_map(tris, E, uE, edgeUeInfo, uE2E);
  const int num_faces = tris.rows();
  D.setConstant(tris.rows(),tris.cols(),false);
  // loop over all unique edges
  for(int ue = 0;ue < uE2E.size(); ue++)
  {
    const bool ue_is_d = is_delaunay(vers,tris,uE2E,float_incircle,ue);
    // Set for all instances
    for(int e = 0;e<uE2E[ue].size();e++)
    {
      D( uE2E[ue][e]%num_faces, uE2E[ue][e]/num_faces) = ue_is_d;
    }
  }
}

template <
  typename DerivedV,
  typename DerivedF,
  typename uE2EType,
  typename InCircle,
  typename ueiType>
IGL_INLINE bool igl::is_delaunay(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  const std::vector<std::vector<uE2EType> > & uE2E,
  const InCircle incircle,
  const ueiType uei)
{
  if(uE2E[uei].size() == 1) return true;
  if(uE2E[uei].size() > 2) return false;
  const int num_faces = tris.rows();
  typedef typename DerivedV::Scalar Scalar;
  const auto& half_edges = uE2E[uei];
  assert((half_edges.size() == 2)  && "uE2E[uei].size() should be 2");
  const size_t f1 = half_edges[0] % num_faces;
  const size_t f2 = half_edges[1] % num_faces;
  const size_t c1 = half_edges[0] / num_faces;
  const size_t c2 = half_edges[1] / num_faces;
  assert(c1 < 3);
  assert(c2 < 3);
  assert(f1 != f2);
  const size_t v1 = tris(f1, (c1+1)%3);
  const size_t v2 = tris(f1, (c1+2)%3);
  const size_t v4 = tris(f1, c1);
  const size_t v3 = tris(f2, c2);
  const Scalar p1[] = {vers(v1, 0), vers(v1, 1)};
  const Scalar p2[] = {vers(v2, 0), vers(v2, 1)};
  const Scalar p3[] = {vers(v3, 0), vers(v3, 1)};
  const Scalar p4[] = {vers(v4, 0), vers(v4, 1)};
  auto orientation = incircle(p1, p2, p4, p3);
  return orientation <= 0;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::is_delaunay<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<bool, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 3, 0, -1, 3> >&);
// generated by autoexplicit.sh
template void igl::is_delaunay<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<bool, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::is_delaunay<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, int, short (*)(double const*, double const*, double const*, double const*), unsigned long>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, short (*)(double const*, double const*, double const*, double const*), unsigned long);
#ifdef WIN32
template bool igl::is_delaunay<class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, int, short(*)(double const *, double const *, double const *, double const *), unsigned __int64>(class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class std::vector<class std::vector<int, class std::allocator<int>>, class std::allocator<class std::vector<int, class std::allocator<int>>>> const &, short(*const)(double const *, double const *, double const *, double const *), unsigned __int64);
#endif
#endif
