#include "pseudonormal_test.h"
#include "barycentric_coordinates.h"
#include "doublearea.h"
#include "project_to_line_segment.h"
#include <cassert>


template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedFN,
  typename DerivedVN,
  typename DerivedEN,
  typename DerivedEMAP,
  typename Derivedq,
  typename Derivedc,
  typename Scalar,
  typename Derivedn>
IGL_INLINE void igl::pseudonormal_test(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  const Eigen::MatrixBase<DerivedFN> & triNorms,
  const Eigen::MatrixBase<DerivedVN> & verNorms,
  const Eigen::MatrixBase<DerivedEN> & edgeNorms,
  const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo,
  const Eigen::MatrixBase<Derivedq> & q,
  const int triIdx0,
  Eigen::PlainObjectBase<Derivedc> & pos0,
  Scalar & sign,
  Eigen::PlainObjectBase<Derivedn> & normal)
{
  using namespace Eigen;
  const auto & qc = q-pos0;
  typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
  RowVector3S b;
  // Using barycentric coorindates to determine whether close to a vertex/edge
  // seems prone to error when dealing with nearly degenerate triangles: Even
  // the barycenter (1/3,1/3,1/3) can be made arbitrarily close to an
  // edge/vertex
  //
  const RowVector3S A = vers.row(tris(triIdx0,0));
  const RowVector3S B = vers.row(tris(triIdx0,1));
  const RowVector3S C = vers.row(tris(triIdx0,2));

  const double area = [&A,&B,&C]()
  {
    Matrix<double,1,1> area;
    doublearea(A,B,C,area);
    return area(0);
  }();
  // These were chosen arbitrarily. In a floating point scenario, I'm not sure
  // the best way to determine if pos0 is on a vertex/edge or in the middle of the
  // face: specifically, I'm worrying about degenerate triangles where
  // barycentric coordinates are error-prone.
  const double MIN_DOUBLE_AREA = 1e-4;
  const double epsilon = 1e-12;
  if(area>MIN_DOUBLE_AREA)
  {
    barycentric_coordinates( pos0,A,B,C,b);
    // Determine which normal to use
    const int type = (b.array()<=epsilon).template cast<int>().sum();
    switch(type)
    {
      case 2:
        // Find vertex
        for(int x = 0;x<3;x++)
        {
          if(b(x)>epsilon)
          {
            normal = verNorms.row(tris(triIdx0,x));
            break;
          }
        }
        break;

      case 1:
        // Find edge
        for(int x = 0;x<3;x++)
        {
          if(b(x)<=epsilon)
          {
            normal = edgeNorms.row(edgeUeInfo(tris.rows()*x+triIdx0));
            break;
          }
        }
        break;
      default:
        assert(false && "all barycentric coords zero.");
      case 0:
        normal = triNorms.row(triIdx0);
        break;
    }
  }
  else
  {
    // Check each vertex
    bool found = false;
    for(int v = 0;v<3 && !found;v++)
    {
      if( (pos0-vers.row(tris(triIdx0,v))).norm() < epsilon)
      {
        found = true;
        normal = verNorms.row(tris(triIdx0,v));
      }
    }

    // Check each edge
    for(int e = 0;e<3 && !found;e++)
    {
      const RowVector3S sign = vers.row(tris(triIdx0,(e+1)%3));
      const RowVector3S d = vers.row(tris(triIdx0,(e+2)%3));
      Matrix<double,1,1> sqr_d_j_x(1,1);
      Matrix<double,1,1> t(1,1);
      project_to_line_segment(pos0,sign,d,t,sqr_d_j_x);
      if(sqrt(sqr_d_j_x(0)) < epsilon)
      {
        normal = edgeNorms.row(edgeUeInfo(tris.rows()*e+triIdx0));
        found = true;
      }
    }

    // Finally just use face
    if(!found) 
      normal = triNorms.row(triIdx0); 
  }
  sign = (qc.dot(normal) >= 0 ? 1. : -1.);
}



// V-EÕ¯∏Ò π”√£∫
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedEN,
  typename DerivedVN,
  typename Derivedq,
  typename Derivedc,
  typename Scalar,
  typename Derivedn>
IGL_INLINE void igl::pseudonormal_test(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & E,
  const Eigen::MatrixBase<DerivedEN> & edgeNorms,
  const Eigen::MatrixBase<DerivedVN> & verNorms,
  const Eigen::MatrixBase<Derivedq> & q,
  const int e,
  Eigen::PlainObjectBase<Derivedc> & pos0,
  Scalar & sign,
  Eigen::PlainObjectBase<Derivedn> & normal)
{
  using namespace Eigen;
  const auto & qc = q-pos0;
  const double len = (vers.row(E(e,1))-vers.row(E(e,0))).norm();
  // barycentric coordinates
  // this .head() nonsense is for "ridiculus" templates instantiations that AABB needs to compile
  Eigen::Matrix<Scalar,1,2>
    b((pos0-vers.row(E(e,1))).norm()/len,(pos0-vers.row(E(e,0))).norm()/len);
    //b((pos0-vers.row(E(e,1)).head(pos0.size())).norm()/len,(pos0-vers.row(E(e,0)).head(pos0.size())).norm()/len);
  // Determine which normal to use
  const double epsilon = 1e-12;
  const int type = (b.array()<=epsilon).template cast<int>().sum();
  switch(type)
  {
    case 1:
      // Find vertex
      for(int x = 0;x<2;x++)
      {
        if(b(x)>epsilon)
        {
          normal = verNorms.row(E(e,x)).head(2);
          break;
        }
      }
      break;
    default:
      assert(false && "all barycentric coords zero.");
    case 0:
      normal = edgeNorms.row(e).head(2);
      break;
  }
  sign = (qc.dot(normal) >= 0 ? 1. : -1.);
}


// This is a bullshit template because AABB annoyingly needs templates for bad combinations of 3D vers with DIM=2 AABB
// _Define_ as a no-op rather than monkeying around with the proper code above
namespace igl
{
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&) {assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&){assert(false);};
  template <> IGL_INLINE void pseudonormal_test(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&){assert(false);};
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
// NEW
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, float, Eigen::Matrix<float, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, float, Eigen::Matrix<float, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, 1, 2, 1, 1, 2>, Eigen::Matrix<float, 1, 2, 1, 1, 2>, float, Eigen::Matrix<float, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, 1, 2, 1, 1, 2>, Eigen::Matrix<float, 1, 2, 1, 1, 2>, float, Eigen::Matrix<float, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
// OLD
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, float, Eigen::Matrix<float, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, Eigen::Matrix<float, 1, -1, 1, 1, -1>, float, Eigen::Matrix<float, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, -1, 1, 1, -1> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, 1, 2, 1, 1, 2>, Eigen::Matrix<float, 1, -1, 1, 1, -1>, float, Eigen::Matrix<float, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, -1, 1, 1, -1> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, -1, 1, 1, -1>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, 1, -1, 1, 1, -1>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 1, -1, false>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 1, -1, false> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<float, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, Eigen::Matrix<float, 1, 3, 1, 1, 3>, float, Eigen::Matrix<float, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&, float&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3>, double, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template void igl::pseudonormal_test<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, 1, 2, 1, 1, 2>, double, Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, double&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&);
#endif
