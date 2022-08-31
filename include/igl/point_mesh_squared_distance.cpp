#include "point_mesh_squared_distance.h"
#include "AABB.h"
#include <cassert>


// point_mesh_squared_distance()——计算点到网格的距离；
template <
  typename DerivedP,
  typename DerivedV,
  typename DerivedEle,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE void igl::point_mesh_squared_distance(
  const Eigen::MatrixBase<DerivedP> & vers0,
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  Eigen::PlainObjectBase<DerivedsqrD> & minSqrDis,
  Eigen::PlainObjectBase<DerivedI> & minDisIdx,
  Eigen::PlainObjectBase<DerivedC> & minDisVers)
{
  using namespace std;
  const size_t dim = vers0.cols();
  assert((dim == 2 || dim == 3) && "vers0.cols() should be 2 or 3");
  assert(vers0.cols() == vers.cols() && "vers0.cols() should equal vers.cols()");

  switch(dim)
  {
    default:
      assert(false && "Unsupported dim");

      // fall-through and pray
    case 3:
    {
      AABB<DerivedV,3> tree;
      tree.init(vers, Ele);
      return tree.squared_distance(vers, Ele, vers0, minSqrDis, minDisIdx, minDisVers);
    }
    case 2:
    {
      AABB<DerivedV,2> tree;
      tree.init(vers,Ele);
      return tree.squared_distance(vers,Ele,vers0,minSqrDis,minDisIdx,minDisVers);
    }
  }
}



// 模板特化：
#ifdef IGL_STATIC_LIBRARY
template void igl::point_mesh_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> &);
template void igl::point_mesh_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::point_mesh_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&);
template void igl::point_mesh_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::point_mesh_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);

#ifdef WIN32
template void igl::point_mesh_squared_distance<class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<__int64, -1, 1, 0, -1, 1>, class Eigen::Matrix<double, -1, 3, 0, -1, 3>>(class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, 1, 0, -1, 1>> &, class Eigen::PlainObjectBase<class Eigen::Matrix<__int64, -1, 1, 0, -1, 1>> &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, 3, 0, -1, 3>> &);
template void igl::point_mesh_squared_distance<class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<double, -1, 1, 0, -1, 1>, class Eigen::Matrix<__int64, -1, 1, 0, -1, 1>, class Eigen::Matrix<double, -1, 3, 0, -1, 3>>(class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1>> const &, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1>> const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, 1, 0, -1, 1>> &, class Eigen::PlainObjectBase<class Eigen::Matrix<__int64, -1, 1, 0, -1, 1>> &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, 3, 0, -1, 3>> &);
#endif
#endif
