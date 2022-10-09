#ifndef IGL_WRITEOBJ_H
#define IGL_WRITEOBJ_H
#include "igl_inline.h"
// History:
//  return type changed from void to bool  Alec 20 Sept 2011

#include <Eigen/Core>
#include <string>
#include <vector>

namespace igl 
{
  // Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   vers  #vers by 3 mesh vertex positions
  //   tris  #tris by 3|4 mesh indices into vers
  //   CN #CN by 3 normal vectors
  //   FN  #tris by 3|4 corner normal indices into CN
  //   TC  #TC by 2|3 texture coordinates
  //   FTC #tris by 3|4 corner texture coord indices into TC
  // Returns true on success, false on error
  //
  // Known issues: Horrifyingly, this does not have the same order of
  // parameters as readOBJ.
  template <
    typename DerivedV, 
    typename DerivedF,
    typename DerivedCN, 
    typename DerivedFN,
    typename DerivedTC, 
    typename DerivedFTC>
  IGL_INLINE bool writeOBJ(
    const std::string str,
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const Eigen::MatrixBase<DerivedCN>& CN,
    const Eigen::MatrixBase<DerivedFN>& FN,
    const Eigen::MatrixBase<DerivedTC>& TC,
    const Eigen::MatrixBase<DerivedFTC>& FTC);
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeOBJ(
    const std::string str,
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris);

  // Write a mesh of mixed tris and quads to an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   vers  #vers by 3 mesh vertex positions
  //   tris  #tris std::vector of std::vector<Index> of size 3 or 4 mesh indices into vers
  // Returns true on success, false on error
  template <typename DerivedV, typename T>
  IGL_INLINE bool writeOBJ(
    const std::string &str,
    const Eigen::MatrixBase<DerivedV>& vers,
    const std::vector<std::vector<T> >& tris);

}

#ifndef IGL_STATIC_LIBRARY
#  include "writeOBJ.cpp"
#endif

#endif
