#include "unique_edge_map.h"
#include "oriented_facets.h"
#include "unique_simplices.h"
#include "cumsum.h"
#include "accumarray.h"
#include <cassert>
#include <algorithm>


// 重载1.1
template <
  typename DerivedF, 
  typename DerivedE, 
  typename DeriveduE, 
  typename DerivedEMAP, 
  typename uE2EType>
IGL_INLINE void igl::unique_edge_map(
  const Eigen::MatrixBase<DerivedF> & tris, 
  Eigen::PlainObjectBase<DerivedE> & edges, 
  Eigen::PlainObjectBase<DeriveduE> & uE, 
  Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo, 
  std::vector<std::vector<uE2EType> > & uE2E)
{
  using namespace Eigen;
  using namespace std;

  unique_edge_map(tris, edges, uE, edgeUeInfo);
  uE2E.resize(uE.rows());

  // This does help a little
  for_each(uE2E.begin(), uE2E.end(), [](vector<uE2EType > & v){v.reserve(2);});
  const size_t edgesCount = edges.rows();
  assert((size_t)edgeUeInfo.size() == edgesCount);
  for(uE2EType e = 0;e<(uE2EType)edgesCount;e++)
        uE2E[edgeUeInfo(e)].push_back(e);
}


// 重载1.
template <
  typename DerivedF, 
  typename DerivedE, 
  typename DeriveduE, 
  typename DerivedEMAP>
IGL_INLINE void igl::unique_edge_map(
  const Eigen::MatrixBase<DerivedF> & tris, 
  Eigen::PlainObjectBase<DerivedE> & edges, 
  Eigen::PlainObjectBase<DeriveduE> & uE, 
  Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo)
{
  using namespace Eigen;
  using namespace std;

  oriented_facets(tris,  edges);                        // 计算网格的半边（有向边 ）
  const size_t edgesCount = edges.rows();

  // This is 2x faster to create than a map from pairs to lists of edges and 5x faster to access
  //  (actually access is probably assympotically faster O(1) vs. O(log m)
  Matrix<typename DerivedEMAP::Scalar, Dynamic, 1> IA;
  unique_simplices(edges, uE, IA, edgeUeInfo);

  assert((size_t)edgeUeInfo.size() == edgesCount);
}


// 重载1.2
template <
  typename DerivedF, 
  typename DerivedE, 
  typename DeriveduE, 
  typename DerivedEMAP, 
  typename DeriveduEC, 
  typename DeriveduEE>
IGL_INLINE void igl::unique_edge_map(
  const Eigen::MatrixBase<DerivedF> & tris, 
  Eigen::PlainObjectBase<DerivedE> & edges, 
  Eigen::PlainObjectBase<DeriveduE> & uE, 
  Eigen::PlainObjectBase<DerivedEMAP> & edgeUeInfo, 
  Eigen::PlainObjectBase<DeriveduEC> & uEC, 
  Eigen::PlainObjectBase<DeriveduEE> & uEE)
{
  // Avoid using uE2E
  igl::unique_edge_map(tris, edges, uE, edgeUeInfo);
  assert(edgeUeInfo.maxCoeff() < uE.rows());

  // counts of each unique edge
  typedef Eigen::Matrix<typename DeriveduEC::Scalar, Eigen::Dynamic, 1> VectorXI;
  VectorXI uEK;
  igl::accumarray(edgeUeInfo, 1, uEK);
  assert(uEK.rows() == uE.rows());

  // base offset in uEE
  igl::cumsum(uEK, 1, true, uEC);
  assert(uEK.rows()+1 == uEC.rows());

  // running inner offset in uEE
  VectorXI uEO = VectorXI::Zero(uE.rows(), 1);

  // flat array of faces incide on each uE
  uEE.resize(edgeUeInfo.rows(), 1);
  for(Eigen::Index e = 0;e<edgeUeInfo.rows();e++)
  {
    const typename DerivedEMAP::Scalar ue = edgeUeInfo(e);
    const typename DeriveduEC::Scalar i = uEC(ue)+ uEO(ue);
    uEE(i) = e;
    uEO(ue)++;
  }

  assert( (uEK.array()==uEO.array()).all() );
}



#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::unique_edge_map<\
            Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  \
            Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  \
            Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  \
            Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  \
            Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  \
            Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(\
        Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  \
        Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  \
        Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  \
        Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  \
        Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  \
        Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);


// generated by autoexplicit.sh
template void igl::unique_edge_map<Eigen::Matrix<unsigned int,  -1,  -1,  1,  -1,  -1>,  Eigen::Matrix<unsigned int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<unsigned int,  -1,  -1,  1,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<unsigned int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);
// generated by autoexplicit.sh
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);
// generated by autoexplicit.sh
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);
// generated by autoexplicit.sh
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  unsigned long>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<unsigned long,  std::allocator<unsigned long> >,  std::allocator<std::vector<unsigned long,  std::allocator<unsigned long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  int>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  std::vector<std::vector<int,  std::allocator<int> >,  std::allocator<std::vector<int,  std::allocator<int> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  unsigned long>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  std::vector<std::vector<unsigned long,  std::allocator<unsigned long> >,  std::allocator<std::vector<unsigned long,  std::allocator<unsigned long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<long,  -1,  1,  0,  -1,  1>,  long>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<long,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<long,  std::allocator<long> >,  std::allocator<std::vector<long,  std::allocator<long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<long,  -1,  1,  0,  -1,  1>,  long>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<long,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<long,  std::allocator<long> >,  std::allocator<std::vector<long,  std::allocator<long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  int>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<int,  std::allocator<int> >,  std::allocator<std::vector<int,  std::allocator<int> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  int>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<int,  std::allocator<int> >,  std::allocator<std::vector<int,  std::allocator<int> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<double,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<long,  -1,  1,  0,  -1,  1>,  long>(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<long,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<long,  std::allocator<long> >,  std::allocator<std::vector<long,  std::allocator<long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  int>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<int,  std::allocator<int> >,  std::allocator<std::vector<int,  std::allocator<int> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  unsigned long>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<unsigned long,  std::allocator<unsigned long> >,  std::allocator<std::vector<unsigned long,  std::allocator<unsigned long> > > >&);
template void igl::unique_edge_map<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  int>(Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  2,  0,  -1,  2> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  std::vector<std::vector<int,  std::allocator<int> >,  std::allocator<std::vector<int,  std::allocator<int> > > >&);

#ifdef WIN32
template void igl::unique_edge_map<class Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  class Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  class Eigen::Matrix<int,  -1,  2,  0,  -1,  2>,  class Eigen::Matrix<__int64,  -1,  1,  0,  -1,  1>,  __int64>(class Eigen::MatrixBase<class Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  2,  0,  -1,  2> > &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  2,  0,  -1,  2> > &,  class Eigen::PlainObjectBase<class Eigen::Matrix<__int64,  -1,  1,  0,  -1,  1> > &,  class std::vector<class std::vector<__int64,  class std::allocator<__int64> >,  class std::allocator<class std::vector<__int64,  class std::allocator<__int64> > > > &);
template void igl::unique_edge_map<class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, 2, 0, -1, 2>, class Eigen::Matrix<int, -1, 2, 0, -1, 2>, class Eigen::Matrix<__int64, -1, 1, 0, -1, 1>, __int64>(class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 2, 0, -1, 2> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 2, 0, -1, 2> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<__int64, -1, 1, 0, -1, 1> > &, class std::vector<class std::vector<__int64, class std::allocator<__int64> >, class std::allocator<class std::vector<__int64, class std::allocator<__int64> > > > &);
template void igl::unique_edge_map<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  class Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  unsigned __int64>(class Eigen::MatrixBase<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>> const &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>> &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>> &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  1,  0,  -1,  1>> &,  class std::vector<class std::vector<unsigned __int64,  class std::allocator<unsigned __int64>>,  class std::allocator<class std::vector<unsigned __int64,  class std::allocator<unsigned __int64>>>> &);
template void igl::unique_edge_map<class Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  class Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  unsigned __int64>(class Eigen::MatrixBase<class Eigen::Matrix<int,  -1,  3,  1,  -1,  3>> const &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>> &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>> &,  class Eigen::PlainObjectBase<class Eigen::Matrix<int,  -1,  1,  0,  -1,  1>> &,  class std::vector<class std::vector<unsigned __int64,  class std::allocator<unsigned __int64>>,  class std::allocator<class std::vector<unsigned __int64,  class std::allocator<unsigned __int64>>>> &);
template void igl::unique_edge_map<class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned __int64>(class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, class std::vector<class std::vector<unsigned __int64, class std::allocator<unsigned __int64> >, class std::allocator<class std::vector<unsigned __int64, class std::allocator<unsigned __int64> > > > &);
#endif

#endif
