#ifndef IGL_READMESH_H
#define IGL_READMESH_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace igl
{
    // readMESH()——读取.mesh四面体体素网格文件；
    /*
       load a tetrahedral volume mesh from a .mesh file
  
       Templates:
         Scalar  type for positions and vectors (will be read as double and cast
           to Scalar)
         Index  type for indices (will be read as int and cast to Index)
       
       Input:
         mesh_file_name  path of .mesh file
       
       Outputs:
         V  double matrix of vertex positions  #V by 3
         T  #T list of tet indices into vertex positions
         F  #F list of face indices into vertex positions
  
       Known bugs: Holes and regions are not supported
    */
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & F);
  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & F);

  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   T  eigen int matrix #T by 4
  //   F  eigen int matrix #F by 3
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F);
  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F);
}

#ifndef IGL_STATIC_LIBRARY
#  include "readMESH.cpp"
#endif

#endif
