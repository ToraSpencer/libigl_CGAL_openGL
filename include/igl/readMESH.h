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
         vers           double matrix of vertex positions  #vers by 3
         T              #T list of tet indices into vertex positions
         tris            #tris list of face indices into vertex positions
  
       Known bugs: Holes and regions are not supported
    */
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    std::vector<std::vector<Scalar > > & vers,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & tris);


  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename Scalar, typename Index>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    std::vector<std::vector<Scalar > > & vers,
    std::vector<std::vector<Index > > & T,
    std::vector<std::vector<Index > > & tris);

  // Input:
  //   mesh_file_name  path of .mesh file
  // Outputs:
  //   vers  eigen double matrix #vers by 3
  //   T  eigen int matrix #T by 4
  //   tris  eigen int matrix #tris by 3
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    const std::string mesh_file_name,
    Eigen::PlainObjectBase<DerivedV>& vers,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& tris);


  // Inputs:
  //   mesh_file  pointer to already opened .mesh file 
  // Outputs:
  //   mesh_file  closed file
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readMESH(
    FILE * mesh_file,
    Eigen::PlainObjectBase<DerivedV>& vers,
    Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& tris);
}

#ifndef IGL_STATIC_LIBRARY
#  include "readMESH.cpp"
#endif

#endif
