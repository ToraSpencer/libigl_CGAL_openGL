#ifndef IGL_SIGNED_DISTANCE_H
#define IGL_SIGNED_DISTANCE_H

// 计算网格符号距离场：
#include "igl_inline.h"
#include "AABB.h"
#include "WindingNumberAABB.h"
#include "fast_winding_number.h"
#include <Eigen/Core>
#include <vector>


namespace igl
{
  enum SignedDistanceType
  {
    // Use fast pseudo-normal test [Bærentzen & Aanæs 2005]
    SIGNED_DISTANCE_TYPE_PSEUDONORMAL         = 0,

    // Use winding number [Jacobson, Kavan Sorking-Hornug 2013]
    SIGNED_DISTANCE_TYPE_WINDING_NUMBER       = 1,
    SIGNED_DISTANCE_TYPE_DEFAULT              = 2,
    SIGNED_DISTANCE_TYPE_UNSIGNED             = 3,

    // Use Fast winding number [Barill, Dickson, Schmidt, Levin, Jacobson 2018]
    SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER  = 4,
    NUM_SIGNED_DISTANCE_TYPE                  = 5
  };


  /*
     Computes signed distance to a mesh
  
       Inputs:
         P  #P by 3 list of query point positions
         vers  #vers by 3 list of vertex positions
         tris  #tris by ss list of triangle indices, ss should be 3 unless sign_type ==
           SIGNED_DISTANCE_TYPE_UNSIGNED
         sign_type  method for computing distance _sign_ SDF
         lower_bound  lower bound of distances needed {std::numeric_limits::min}
         upper_bound  lower bound of distances needed {std::numeric_limits::max}
  
       Outputs:
         SDF  #P list of smallest signed distances
         I  #P list of facet indices corresponding to smallest distances
         C  #P by 3 list of closest points
         N  #P by 3 list of closest normals (only set if sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
  
       Known bugs: This only computes distances to triangles. So unreferenced
       vertices and degenerate triangles are ignored.
  */
  template <
    typename DerivedP,
    typename DerivedV,
    typename DerivedF,
    typename DerivedS,
    typename DerivedI,
    typename DerivedC,
    typename DerivedN>
  IGL_INLINE void signed_distance(
    const Eigen::MatrixBase<DerivedP> & P,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const SignedDistanceType sign_type,
    const typename DerivedV::Scalar lower_bound,
    const typename DerivedV::Scalar upper_bound,
    Eigen::PlainObjectBase<DerivedS> & SDF,
    Eigen::PlainObjectBase<DerivedI> & I,
    Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Computes signed distance to a mesh, with default bounds
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   vers  #vers by 3 list of vertex positions
  //   tris  #tris by ss list of triangle indices, ss should be 3 unless sign_type ==
  //     SIGNED_DISTANCE_TYPE_UNSIGNED
  //   sign_type  method for computing distance _sign_ SDF
  //   lower_bound  lower bound of distances needed {std::numeric_limits::min}
  //   upper_bound  lower bound of distances needed {std::numeric_limits::max}
  // Outputs:
  //   SDF  #P list of smallest signed distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //   N  #P by 3 list of closest normals (only set if
  //     sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
  template <
    typename DerivedP,
    typename DerivedV,
    typename DerivedF,
    typename DerivedS,
    typename DerivedI,
    typename DerivedC,
    typename DerivedN>
  IGL_INLINE void signed_distance(
    const Eigen::MatrixBase<DerivedP> & P,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const SignedDistanceType sign_type,
    Eigen::PlainObjectBase<DerivedS> & SDF,
    Eigen::PlainObjectBase<DerivedI> & I,
    Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Computes signed distance to mesh using pseudonormal with precomputed AABB tree and edge/vertice normals
  //
  // Inputs:
  //   tree  AABB acceleration tree (see AABB.h)
  //   tris  #tris by 3 list of triangle indices
  //   FN  #tris by 3 list of triangle normals 
  //   VN  #vers by 3 list of vertex normals (ANGLE WEIGHTING)
  //   EN  #E by 3 list of edge normals (UNIFORM WEIGHTING)
  //   edgeUeInfo  #tris*3 mapping edges in tris to E
  //   q  Query point
  // Returns signed distance to mesh
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedFN,
    typename DerivedVN,
    typename DerivedEN,
    typename DerivedEMAP,
    typename Derivedq>
  IGL_INLINE typename DerivedV::Scalar signed_distance_pseudonormal(
    const AABB<DerivedV,3> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedFN> & FN,
    const Eigen::MatrixBase<DerivedVN> & VN,
    const Eigen::MatrixBase<DerivedEN> & EN,
    const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo,
    const Eigen::MatrixBase<Derivedq> & q);


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
  IGL_INLINE void signed_distance_pseudonormal(
    const Eigen::MatrixBase<DerivedP> & P,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const AABB<DerivedV,3> & tree,
    const Eigen::MatrixBase<DerivedFN> & FN,
    const Eigen::MatrixBase<DerivedVN> & VN,
    const Eigen::MatrixBase<DerivedEN> & EN,
    const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo,
    Eigen::PlainObjectBase<DerivedS> & SDF,
    Eigen::PlainObjectBase<DerivedI> & I,
    Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Outputs:
  //   s  sign
  //   sqrd  squared distance
  //   i  closest primitive
  //   c  closest point
  //   n  normal
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
  IGL_INLINE void signed_distance_pseudonormal(
    const AABB<DerivedV,3> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedFN> & FN,
    const Eigen::MatrixBase<DerivedVN> & VN,
    const Eigen::MatrixBase<DerivedEN> & EN,
    const Eigen::MatrixBase<DerivedEMAP> & edgeUeInfo,
    const Eigen::MatrixBase<Derivedq> & q,
    Scalar & s,
    Scalar & sqrd,
    int & i,
    Eigen::PlainObjectBase<Derivedc> & c,
    Eigen::PlainObjectBase<Derivedn> & n);


  template <
    typename DerivedV,
    typename DerivedE,
    typename DerivedEN,
    typename DerivedVN,
    typename Derivedq,
    typename Scalar,
    typename Derivedc,
    typename Derivedn>
  IGL_INLINE void signed_distance_pseudonormal(
    const AABB<DerivedV,2> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedE> & E,
    const Eigen::MatrixBase<DerivedEN> & EN,
    const Eigen::MatrixBase<DerivedVN> & VN,
    const Eigen::MatrixBase<Derivedq> & q,
    Scalar & s,
    Scalar & sqrd,
    int & i,
    Eigen::PlainObjectBase<Derivedc> & c,
    Eigen::PlainObjectBase<Derivedn> & n);


  // Inputs:
  //   tree  AABB acceleration tree (see cgal/point_mesh_squared_distance.h)
  //   hier  Winding number evaluation hierarchy
  //   q  Query point
  // Returns signed distance to mesh
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedq>
  IGL_INLINE typename DerivedV::Scalar signed_distance_winding_number(
    const AABB<DerivedV,3> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const igl::WindingNumberAABB<Derivedq,DerivedV,DerivedF> & hier,
    const Eigen::MatrixBase<Derivedq> & q);


  // Outputs:
  //   s  sign
  //   sqrd  squared distance
  //   pp  closest point and primitve
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedq,
    typename Scalar,
    typename Derivedc>
  IGL_INLINE void signed_distance_winding_number(
    const AABB<DerivedV,3> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const igl::WindingNumberAABB<Derivedq,DerivedV,DerivedF> & hier,
    const Eigen::MatrixBase<Derivedq> & q,
    Scalar & s,
    Scalar & sqrd,
    int & i,
    Eigen::PlainObjectBase<Derivedc> & c);


  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedq,
    typename Scalar,
    typename Derivedc>
  IGL_INLINE void signed_distance_winding_number(
    const AABB<DerivedV,2> & tree,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<Derivedq> & q,
    Scalar & s,
    Scalar & sqrd,
    int & i,
    Eigen::PlainObjectBase<Derivedc> & c);


  // Calculates signed distance at query points P, using fast winding number
  //   for sign.
  //
  // Usage:
  //     VectorXd SDF;  
  //     VectorXd vers, P; //where vers is mesh vertices, P are query points
  //     VectorXi tris;  
  //     igl::FastWindingNumberBVH fwn_bvh;
  //     igl::fast_winding_number(vers.cast<float>(), tris, 2, fwn_bvh);
  //     igl::signed_distance_fast_winding_number(P,vers,tris,tree,fwn_bvh,SDF);
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   vers  #vers by 3 list of triangle indices
  //   tris  #tris by 3 list of triangle normals 
  //   tree  AABB acceleration tree (see AABB.h)
  //   bvh fast winding precomputation (see Fast_Winding_Number.h)   
  // Outputs:
  //   SDF  #P list of signed distances of each point in P
  template <
    typename DerivedP,
    typename DerivedV,
    typename DerivedF,
    typename DerivedS>
  IGL_INLINE void signed_distance_fast_winding_number(
    const Eigen::MatrixBase<DerivedP> & P,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const AABB<DerivedV,3> & tree,
    const igl::FastWindingNumberBVH & fwn_bvh,
    Eigen::PlainObjectBase<DerivedS> & SDF
  );

  // Calculates signed distance at query point q, using fast winding number
  //   for sign.
  //
  // Inputs:
  //   tree  AABB acceleration tree (see AABB.h)
  //   vers  #vers by 3 list of triangle indices
  //   tris  #tris by 3 list of triangle normals 
  //   bvh fast winding precomputation (see Fast_Winding_Number.h)   
  //   q  1 by 3 list of query point positions
  // Outputs:
  //   SDF  #P list of signed distances of each point in P
  template <
    typename Derivedq,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE typename DerivedV::Scalar signed_distance_fast_winding_number(
    const Eigen::MatrixBase<Derivedq> & q,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const AABB<DerivedV,3> & tree,
    const igl::FastWindingNumberBVH & fwn_bvh
  );
}

#ifndef IGL_STATIC_LIBRARY
#  include "signed_distance.cpp"
#endif

#endif
