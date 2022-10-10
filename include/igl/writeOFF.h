#ifndef IGL_WRITEOFF_H
#define IGL_WRITEOFF_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
  //Export geometry and colors-by-vertex
  // Export a mesh from an ascii OFF file, filling in vertex positions.
  // Only triangle meshes are supported
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .off output file
  //   vers  #vers by 3 mesh vertex positions
  //   tris  #tris by 3 mesh indices into vers
  //   C  double matrix of rgb values per vertex #vers by 3
  // Outputs:
  // Returns true on success, false on errors
  template <typename DerivedV, typename DerivedF, typename DerivedC>
  IGL_INLINE bool writeOFF(
    const std::string str,
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const Eigen::MatrixBase<DerivedC>& C);

  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeOFF(
    const std::string str,
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeOFF.cpp"
#endif

#endif
