#ifndef IGL_WRITESTL_H
#define IGL_WRITESTL_H
#include "igl_inline.h"
#include <igl/FileEncoding.h>

#ifndef IGL_NO_EIGEN
#  include <Eigen/Core>
#endif
#include <string>
#include <vector>

namespace igl
{
  // Write a mesh to an stl file.
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  // Inputs:
  //   filename path to .obj file
  //   vers  double matrix of vertex positions  #tris*3 by 3
  //   tris  index matrix of triangle indices #tris by 3
  //   N  double matrix of vertex positions  #tris by 3
  //   encoding enum to set file encoding (ascii by default)
  // Returns true on success, false on errors
  //
  template <typename DerivedV, typename DerivedF, typename DerivedN>
  IGL_INLINE bool writeSTL(
    const std::string & filename,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    const Eigen::MatrixBase<DerivedN> & N,
    FileEncoding encoding=FileEncoding::Ascii);
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeSTL(
    const std::string & filename,
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedF> & tris,
    FileEncoding encoding=FileEncoding::Ascii);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeSTL.cpp"
#endif

#endif
