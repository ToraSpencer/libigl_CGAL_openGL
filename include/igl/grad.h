#ifndef IGL_GRAD_H
#define IGL_GRAD_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl {
  // GRAD
  // G = grad(vers,tris)
  //
  // Compute the numerical gradient operator
  //
  // Inputs:
  //   vers          #vertices by 3 list of mesh vertex positions
  //   tris          #faces by 3 list of mesh face indices [or a #faces by 4 list of tetrahedral indices]
  //   uniform    boolean (default false) - Use a uniform mesh instead of the vertices vers
  // Outputs:
  //   G  #faces*dim by #vers Gradient operator
  //

  // Gradient of a scalar function defined on piecewise linear elements (mesh)
  // is constant on each triangle [tetrahedron] i,j,k:
  // grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  // where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  // i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of
  // 90 degrees
  //
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void grad(
    const Eigen::MatrixBase<DerivedV>&vers,
    const Eigen::MatrixBase<DerivedF>&tris,
    Eigen::SparseMatrix<typename DerivedV::Scalar> &G,
    bool uniform = false);
}
#ifndef IGL_STATIC_LIBRARY
#  include "grad.cpp"
#endif

#endif
