#include "signed_distance.h"
#include "get_seconds.h"
#include "per_edge_normals.h"
#include "parallel_for.h"
#include "per_face_normals.h"
#include "per_vertex_normals.h"
#include "point_mesh_squared_distance.h"
#include "pseudonormal_test.h"
#include "fast_winding_number.h"


template <
  typename DerivedP, 
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedS, 
  typename DerivedI, 
  typename DerivedC, 
  typename DerivedN>
IGL_INLINE void igl::signed_distance(
  const Eigen::MatrixBase<DerivedP> & gridCenters, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const SignedDistanceType sign_type, 
  const typename DerivedV::Scalar lower_bound, 
  const typename DerivedV::Scalar upper_bound, 
  Eigen::PlainObjectBase<DerivedS> & SDFvalues, 
  Eigen::PlainObjectBase<DerivedI> & I, 
  Eigen::PlainObjectBase<DerivedC> & C, 
  Eigen::PlainObjectBase<DerivedN> & N)
{
  using namespace Eigen;
  using namespace std;

  const int dim = vers.cols();

  assert((vers.cols() == 3||vers.cols() == 2) && "vers should have 3d or 2d positions");
  assert((gridCenters.cols() == 3||gridCenters.cols() == 2) && "gridCenters should have 3d or 2d positions");
  assert(vers.cols() == gridCenters.cols() && "vers should have same dimension as gridCenters");

  if (sign_type == SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    assert(vers.cols() == 3 && "vers should be 3D for fast winding number");
 
  // Only unsigned distance is supported for non-triangles
  if(sign_type != SIGNED_DISTANCE_TYPE_UNSIGNED) 
    assert(tris.cols() == dim && "tris should have co-dimension 0 simplices"); 

  typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 3> RowVector3S;

  // 1. Prepare distance computation
  AABB<DerivedV, 3> tree3;
  AABB<DerivedV, 2> tree2;
  switch(dim)
  {
    default:
    case 3:
      tree3.init(vers, tris);
      break;
    case 2:
      tree2.init(vers, tris);
      break;
  }

  // 2. 计算所需要的中间数据——面片法线、顶点法线、缠绕数等：
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> FN, VN, EN;
  Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, 2> E;
  Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, 1> edgeUeInfo;
  WindingNumberAABB<RowVector3S, DerivedV, DerivedF> hier3;
  igl::FastWindingNumberBVH fwn_bvh;
  Eigen::VectorXf W;

  switch(sign_type)
  {
    default:
      assert(false && "Unknown SignedDistanceType");

    case SIGNED_DISTANCE_TYPE_UNSIGNED:
      // do nothing
      break;

    case SIGNED_DISTANCE_TYPE_DEFAULT:

    case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
      switch(dim)
      {
        default:
        case 3:
          hier3.set_mesh(vers, tris);
          hier3.grow();
          break;
        case 2:
          // no precomp,  no hierarchy
          break;
      }
      break;

    case SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER:
      // assert above ensures dim == 3
      igl::fast_winding_number(vers.template cast<float>().eval(),  tris,  2,  fwn_bvh);

    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      switch(dim)
      {
        default:
        case 3:
          // "Signed Distance Computation Using the Angle Weighted Pseudonormal" [Bærentzen & Aanæs 2005]
          per_face_normals(vers, tris, FN);
          per_vertex_normals(vers, tris, PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN, VN);
          per_edge_normals(vers, tris, PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN, EN, E, edgeUeInfo);
          break;
        case 2:
          FN.resize(tris.rows(), 2);
          VN = DerivedV::Zero(vers.rows(), 2);
          for(int e = 0;e<tris.rows();e++)
          {
            // rotate edge vector
            FN(e, 0) =  (vers(tris(e, 1), 1)-vers(tris(e, 0), 1));
            FN(e, 1) = -(vers(tris(e, 1), 0)-vers(tris(e, 0), 0));
            FN.row(e).normalize();
            // add to vertex normal
            VN.row(tris(e, 1)) += FN.row(e);
            VN.row(tris(e, 0)) += FN.row(e);
          }
          // normalize to average
          VN.rowwise().normalize();
          break;
      }
      N.resize(gridCenters.rows(), dim);
      break;
  }

  // 3. convert to bounds on (unsiged) squared distances
  typedef typename DerivedV::Scalar Scalar; 
  const Scalar max_abs = std::max(std::abs(lower_bound), std::abs(upper_bound));
  const Scalar up_sqr_d = std::pow(max_abs, 2.0);
  const Scalar low_sqr_d = std::pow(std::max(max_abs-(upper_bound-lower_bound), (Scalar)0.0), 2.0);
  SDFvalues.resize(gridCenters.rows(), 1);
  I.resize(gridCenters.rows(), 1);
  C.resize(gridCenters.rows(), dim);

  // 4. 计算符号距离的并行for循环：
  parallel_for(gridCenters.rows(), [&](const int p)
  {
    RowVector3S ver3D;
    Eigen::Matrix<typename DerivedV::Scalar, 1, 2>  q2;
    switch(gridCenters.cols())
    {
      default:
      case 3:
        ver3D.head(gridCenters.row(p).size()) = gridCenters.row(p);
        break;
      case 2:
        q2 = gridCenters.row(p).head(2);
        break;
    }

    typename DerivedV::Scalar sign=1, sqrd=0;
    Eigen::Matrix<typename DerivedV::Scalar, 1, Eigen::Dynamic>  c;
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> c3;
    Eigen::Matrix<typename DerivedV::Scalar, 1, 2>  c2;
    int triIdx0 = -1;                   // 当前顶点最近的三角片的索引；         

    // pf1. 计算当前顶点和？？？的平方距离：
    sqrd = dim==3 ? tree3.squared_distance(vers, tris, ver3D, low_sqr_d, up_sqr_d, triIdx0, c3):
      tree2.squared_distance(vers, tris, q2, low_sqr_d, up_sqr_d, triIdx0, c2);

    if(sqrd >= up_sqr_d || sqrd < low_sqr_d)
    {
      // Out of bounds gets a nan (nans on grids can be flood filled later using igl::flood_fill)
      SDFvalues(p) = std::numeric_limits<double>::quiet_NaN();
      I(p) = tris.rows()+1;
      C.row(p).setConstant(0);
    }
    else
    {
      // Determine sign
      switch(sign_type)
      {
        default:
          assert(false && "Unknown SignedDistanceType");

        case SIGNED_DISTANCE_TYPE_UNSIGNED:
          break;

        case SIGNED_DISTANCE_TYPE_DEFAULT:

        case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
        {
          Scalar w = 0;
          if(dim == 3)
            sign = 1.-2.*hier3.winding_number(ver3D.transpose());
          else
          {
            assert(!vers.derived().IsRowMajor);
            assert(!tris.derived().IsRowMajor);
            sign = 1.-2.*winding_number(vers, tris, q2);
          }
          break;
        }

        case SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER:
        {
          //assert above ensured 3D
          Scalar w = fast_winding_number(fwn_bvh,  2,  ver3D.template cast<float>().eval());         
          sign = 1.-2.*std::abs(w);  
          break;
        }

        case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
        {
          RowVector3S n3;
          Eigen::Matrix<typename DerivedV::Scalar, 1, 2>  n2;
          dim==3 ? pseudonormal_test(vers, tris, FN, VN, EN, edgeUeInfo, ver3D, triIdx0, c3, sign, n3) :             
            pseudonormal_test(vers, tris, FN, VN, q2, triIdx0, c2, sign, n2);

          Eigen::Matrix<typename DerivedN::Scalar, 1, Eigen::Dynamic>  n;
          (dim==3 ? n = n3.template cast<typename DerivedN::Scalar>() : n = n2.template cast<typename DerivedN::Scalar>());
          N.row(p) = n.template cast<typename DerivedN::Scalar>();
          break;
        }
      }
      I(p) = triIdx0;
      SDFvalues(p) = sign * sqrt(sqrd);
      C.row(p) = (dim==3 ? c=c3 : c=c2).template cast<typename DerivedC::Scalar>();
    }
  }
  , 10000);

}


template <
  typename DerivedP, 
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedS, 
  typename DerivedI, 
  typename DerivedC, 
  typename DerivedN>
IGL_INLINE void igl::signed_distance(
  const Eigen::MatrixBase<DerivedP> & gridCenters, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const SignedDistanceType sign_type, 
  Eigen::PlainObjectBase<DerivedS> & SDFvalues, 
  Eigen::PlainObjectBase<DerivedI> & I, 
  Eigen::PlainObjectBase<DerivedC> & C, 
  Eigen::PlainObjectBase<DerivedN> & N)
{
  typedef typename DerivedV::Scalar Scalar;
  Scalar lower = std::numeric_limits<Scalar>::min();
  Scalar upper = std::numeric_limits<Scalar>::max();
  return signed_distance(gridCenters, vers, tris, sign_type, lower, upper, SDFvalues, I, C, N);
}


template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedFN, 
  typename DerivedVN, 
  typename DerivedEN, 
  typename DerivedEMAP, 
  typename Derivedq>
IGL_INLINE typename DerivedV::Scalar igl::signed_distance_pseudonormal(
  const AABB<DerivedV, 3> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const Eigen::MatrixBase<DerivedFN> & FN, 
  const Eigen::MatrixBase<DerivedVN> & VN, 
  const Eigen::MatrixBase<DerivedEN> & EN, 
  const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo, 
  const Eigen::MatrixBase<Derivedq> & q)
{
  typename DerivedV::Scalar sign, sqrd;
  Eigen::Matrix<typename DerivedV::Scalar, 1, 3> n, c;
  int triIdx0 = -1;
  signed_distance_pseudonormal(tree, vers, tris, FN, VN, EN, edgeUeInfo, q, sign, sqrd, triIdx0, c, n);
  return sign*sqrt(sqrd);
}


template <
  typename DerivedP, 
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedFN, 
  typename DerivedVN, 
  typename DerivedEN, 
  typename DerivedEMAP, 
  typename DerivedS, 
  typename DerivedI, 
  typename DerivedC, 
  typename DerivedN>
IGL_INLINE void igl::signed_distance_pseudonormal(
  const Eigen::MatrixBase<DerivedP> & gridCenters, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const AABB<DerivedV, 3> & tree, 
  const Eigen::MatrixBase<DerivedFN> & FN, 
  const Eigen::MatrixBase<DerivedVN> & VN, 
  const Eigen::MatrixBase<DerivedEN> & EN, 
  const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo, 
  Eigen::PlainObjectBase<DerivedS> & SDFvalues, 
  Eigen::PlainObjectBase<DerivedI> & I, 
  Eigen::PlainObjectBase<DerivedC> & C, 
  Eigen::PlainObjectBase<DerivedN> & N)
{
  using namespace Eigen;
  const size_t np = gridCenters.rows();
  SDFvalues.resize(np, 1);
  I.resize(np, 1);
  N.resize(np, 3);
  C.resize(np, 3);
  typedef typename AABB<DerivedV, 3>::RowVectorDIMS RowVector3S;
# pragma omp parallel for if(np>1000)
  for(std::ptrdiff_t p = 0;p<np;p++)
  {
    typename DerivedV::Scalar sign, sqrd;
    RowVector3S n, c;
    int triIdx0 = -1;
    RowVector3S q = gridCenters.row(p);
    signed_distance_pseudonormal(tree, vers, tris, FN, VN, EN, edgeUeInfo, q, sign, sqrd, triIdx0, c, n);
    SDFvalues(p) = sign*sqrt(sqrd);
    I(p) = triIdx0;
    N.row(p) = n;
    C.row(p) = c;
  }
}


template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedFN, 
  typename DerivedVN, 
  typename DerivedEN, 
  typename DerivedEMAP, 
  typename Derivedq, 
  typename Scalar, 
  typename Derivedc, 
  typename Derivedn>
IGL_INLINE void igl::signed_distance_pseudonormal(
  const AABB<DerivedV, 3> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const Eigen::MatrixBase<DerivedFN> & FN, 
  const Eigen::MatrixBase<DerivedVN> & VN, 
  const Eigen::MatrixBase<DerivedEN> & EN, 
  const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo, 
  const Eigen::MatrixBase<Derivedq> & q, 
  Scalar & sign, 
  Scalar & sqrd, 
  int & triIdx0, 
  Eigen::PlainObjectBase<Derivedc> & c, 
  Eigen::PlainObjectBase<Derivedn> & n)
{
  using namespace Eigen;
  using namespace std;
  //typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 3> RowVector3S;
  // Alec: Why was this constructor around q necessary?
  //sqrd = tree.squared_distance(vers, tris, RowVector3S(q), triIdx0, (RowVector3S&)c);
  // Alec: Why was this constructor around c necessary?
  //sqrd = tree.squared_distance(vers, tris, q, triIdx0, (RowVector3S&)c);
  sqrd = tree.squared_distance(vers, tris, q, triIdx0, c);
  pseudonormal_test(vers, tris, FN, VN, EN, edgeUeInfo, q, triIdx0, c, sign, n);
}


template <
  typename DerivedV, 
  typename DerivedE, 
  typename DerivedEN, 
  typename DerivedVN, 
  typename Derivedq, 
  typename Scalar, 
  typename Derivedc, 
  typename Derivedn>
IGL_INLINE void igl::signed_distance_pseudonormal(
  const AABB<DerivedV, 2> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedE> & E, 
  const Eigen::MatrixBase<DerivedEN> & EN, 
  const Eigen::MatrixBase<DerivedVN> & VN, 
  const Eigen::MatrixBase<Derivedq> & q, 
  Scalar & sign, 
  Scalar & sqrd, 
  int & triIdx0, 
  Eigen::PlainObjectBase<Derivedc> & c, 
  Eigen::PlainObjectBase<Derivedn> & n)
{
  using namespace Eigen;
  using namespace std;
  typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 2> RowVector2S;
  sqrd = tree.squared_distance(vers, E, RowVector2S(q), triIdx0, (RowVector2S&)c);
  pseudonormal_test(vers, E, EN, VN, q, triIdx0, c, sign, n);
}

template <
  typename DerivedV, 
  typename DerivedF, 
  typename Derivedq>
IGL_INLINE typename DerivedV::Scalar igl::signed_distance_winding_number(
  const AABB<DerivedV, 3> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const igl::WindingNumberAABB<Derivedq, DerivedV, DerivedF> & hier, 
  const Eigen::MatrixBase<Derivedq> & q)
{
  typedef typename DerivedV::Scalar Scalar;
  Scalar sign, sqrd;
  Eigen::Matrix<Scalar, 1, 3> c;
  int triIdx0=-1;
  signed_distance_winding_number(tree, vers, tris, hier, q, sign, sqrd, triIdx0, c);
  return sign*sqrt(sqrd);
}


template <
  typename DerivedV, 
  typename DerivedF, 
  typename Derivedq, 
  typename Scalar, 
  typename Derivedc>
IGL_INLINE void igl::signed_distance_winding_number(
  const AABB<DerivedV, 3> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const igl::WindingNumberAABB<Derivedq, DerivedV, DerivedF> & hier, 
  const Eigen::MatrixBase<Derivedq> & q, 
  Scalar & sign, 
  Scalar & sqrd, 
  int & triIdx0, 
  Eigen::PlainObjectBase<Derivedc> & c)
{
  using namespace Eigen;
  using namespace std;
  typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 3> RowVector3S;
  sqrd = tree.squared_distance(vers, tris, RowVector3S(q), triIdx0, (RowVector3S&)c);
  const Scalar w = hier.winding_number(q.transpose());
  sign = 1.-2.*w;
}

template <
  typename DerivedV, 
  typename DerivedF, 
  typename Derivedq, 
  typename Scalar, 
  typename Derivedc>
IGL_INLINE void igl::signed_distance_winding_number(
  const AABB<DerivedV, 2> & tree, 
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  const Eigen::MatrixBase<Derivedq> & q, 
  Scalar & sign, 
  Scalar & sqrd, 
  int & triIdx0, 
  Eigen::PlainObjectBase<Derivedc> & c)
{
  using namespace Eigen;
  using namespace std;
  typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 2> RowVector2S;
  sqrd = tree.squared_distance(vers, tris, RowVector2S(q), triIdx0, (RowVector2S&)c);
  // TODO: using .data() like this is very dangerous... This is assuming
  // colmajor order
  assert(!vers.derived().IsRowMajor);
  assert(!tris.derived().IsRowMajor);
  sign = 1.-2.*winding_number(vers, tris, q);
}

//Multi point by parrallel for on single point
template <
  typename DerivedP, 
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedS>
IGL_INLINE void igl::signed_distance_fast_winding_number(
    const Eigen::MatrixBase<DerivedP> & gridCenters, 
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    const AABB<DerivedV, 3> & tree, 
    const igl::FastWindingNumberBVH & fwn_bvh, 
    Eigen::PlainObjectBase<DerivedS> & SDFvalues)
  {
    typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 3> RowVector3S;
    SDFvalues.resize(gridCenters.rows(), 1);
    int min_parallel = 10000; 
    parallel_for(gridCenters.rows(),  [&](const int p)
    {
      RowVector3S q;
      q.head(gridCenters.row(p).size()) = gridCenters.row(p);
      // get sdf for single point,  update result matrix
      SDFvalues(p) = signed_distance_fast_winding_number(q,  vers,  tris,  tree, fwn_bvh);
    }
    , min_parallel);  
  }

//Single Point
template <
  typename Derivedq, 
  typename DerivedV, 
  typename DerivedF>
IGL_INLINE typename DerivedV::Scalar igl::signed_distance_fast_winding_number(
    const Eigen::MatrixBase<Derivedq> & q, 
    const Eigen::MatrixBase<DerivedV> & vers, 
    const Eigen::MatrixBase<DerivedF> & tris, 
    const AABB<DerivedV, 3> & tree, 
    const igl::FastWindingNumberBVH & fwn_bvh)
  {
    typedef typename DerivedV::Scalar Scalar;
    Scalar sign, sqrd;
    Eigen::Matrix<Scalar, 1, 3> c;
    int triIdx0 = -1;
    sqrd = tree.squared_distance(vers, tris, q, triIdx0, c);
    Scalar w = fast_winding_number(fwn_bvh, 2, q.template cast<float>());
    //0.5 is on surface
    return sqrt(sqrd)*(1.-2.*std::abs(w));
  }



#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::signed_distance<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<double,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  igl::SignedDistanceType,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>::Scalar,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>::Scalar,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&);
// generated by autoexplicit.sh
template void igl::signed_distance<Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  igl::SignedDistanceType,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>::Scalar,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>::Scalar,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&);
// generated by autoexplicit.sh
template void igl::signed_distance<Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<float,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  igl::SignedDistanceType,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>::Scalar,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>::Scalar,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&);
template void igl::signed_distance<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<float,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  igl::SignedDistanceType,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&);
template void igl::signed_distance_pseudonormal<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  1,  3,  1,  1,  3>,  double,  Eigen::Matrix<double,  1,  3,  1,  1,  3>,  Eigen::Matrix<double,  1,  3,  1,  1,  3> >(igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> > const&,  double&,  double&,  int&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> >&);
template void igl::signed_distance<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::SignedDistanceType,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&);
template Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>::Scalar igl::signed_distance_pseudonormal<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  1,  3,  1,  1,  3> >(igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> > const&);
template void igl::signed_distance_pseudonormal<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&);
template void igl::signed_distance<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::SignedDistanceType,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&);
template void igl::signed_distance_winding_number<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  1,  3,  1,  1,  3>,  double,  Eigen::Matrix<double,  1,  3,  1,  1,  3> >(igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::WindingNumberAABB<Eigen::Matrix<double,  1,  3,  1,  1,  3>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> > const&,  double&,  double&,  int&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  1,  3,  1,  1,  3> >&);
template Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>::Scalar igl::signed_distance_winding_number<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  3,  1,  0,  3,  1> >(igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::WindingNumberAABB<Eigen::Matrix<double,  3,  1,  0,  3,  1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  3,  1,  0,  3,  1> > const&);
template void igl::signed_distance_fast_winding_number<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  igl::AABB<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  3> const&,  igl::FastWindingNumberBVH const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  1,  0,  -1,  1> >&);
#endif
