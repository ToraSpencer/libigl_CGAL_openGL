#ifndef IGL_PER_FACE_NORMALS_H
#define IGL_PER_FACE_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>



namespace igl
{
  // Compute face normals via vertex position list, face list
  // Inputs:
  //   vers  #vers by 3 eigen Matrix of mesh vertex 3D positions
  //   tris  #tris by 3 eigen Matrix of face (triangle) indices
  //   Z  3 vector normal given to faces with degenerate normal.
  // Output:
  //   N  #tris by 3 eigen Matrix of mesh face (triangle) 3D normals
  //
  // Example:
  //   // Give degenerate faces (1/3,1/3,1/3)^0.5
  //   per_face_normals(vers,tris,Vector3d(1,1,1).normalized(),N);
  template <typename DerivedV, typename DerivedF, typename DerivedZ, typename DerivedN>
  IGL_INLINE void per_face_normals(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const Eigen::MatrixBase<DerivedZ> & Z,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Wrapper with Z = (0,0,0). Note that this means that row norms will be zero
  // (i.e. not 1) for degenerate normals.
  template <typename DerivedV, typename DerivedF, typename DerivedN>
  IGL_INLINE void per_face_normals(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Special version where order of face indices is guaranteed not to effect
  // output.
  template <typename DerivedV, typename DerivedF, typename DerivedN>
  IGL_INLINE void per_face_normals_stable(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    Eigen::PlainObjectBase<DerivedN> & N);


  // Inputs:
  //   vers  #vers by 3 list of mesh vertex positions
  //   I  #I vectorized list of polygon corner indices into rows of some matrix vers
  //   C  #polygons+1 list of cumulative polygon sizes so that C(i+1)-C(i) = size of
  //     the ith polygon, and so I(C(i)) through I(C(i+1)-1) are the indices of
  //     the ith polygon
  //   corner_threshold  threshold in degrees on sharp angles
  // Outputs:
  //   N  #polygons by 3 list of per face normals
  //   VV  #I+#polygons by 3 list of auxiliary triangle mesh vertex positions
  //   FF  #I by 3 list of triangle indices into rows of VV
  //   J  #I list of indices into original polygons
  template <
    typename DerivedV,
    typename DerivedI,
    typename DerivedC,
    typename DerivedN,
    typename DerivedVV,
    typename DerivedFF,
    typename DerivedJ>
  IGL_INLINE void per_face_normals(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedI> & I,
    const Eigen::MatrixBase<DerivedC> & C,
    Eigen::PlainObjectBase<DerivedN> & N,
    Eigen::PlainObjectBase<DerivedVV> & VV,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedJ> & J);
}

#ifndef IGL_STATIC_LIBRARY
#  include "per_face_normals.cpp"
#endif

#endif
