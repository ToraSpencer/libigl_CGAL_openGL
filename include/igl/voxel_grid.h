 
#ifndef IGL_VOXEL_GRID_H
#define IGL_VOXEL_GRID_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace igl
{

    // 生成三维轴向包围盒box的栅格
  /*
    Construct the cell center positions of a regular voxel grid(lattice) made  of perfectly square voxels.
   
       Inputs:
         box        bounding box to enclose by grid
         s          number of cell centers on largest gridCounts (including 2 * pad_count)
                    跨度最大的那个维度(xyz中的一个)的栅格数；
         pad_count       number of cells beyond box

       Outputs:
         gridCenters  gridCounts(0)*gridCounts(1)*gridCounts(2) by 3 list of cell center positions
         gridCounts  三个维度上栅格的数目；
    */
  template <
    typename Scalar,
    typename DerivedGV,
    typename Derivedside>
  IGL_INLINE void voxel_grid(
    const Eigen::AlignedBox<Scalar,3> & box, 
    const int s,
    const int pad_count,
    Eigen::PlainObjectBase<DerivedGV> & gridCenters,
    Eigen::PlainObjectBase<Derivedside> & gridCounts);


  // 生成输入点云vers的包围盒的栅格；
    /*
 
       Inputs:
         vers        物体点云
         offset     包围盒边缘到物体最小距离；
         s          number of cell centers on largest gridCounts (including 2 * pad_count)
                    跨度最大的那个维度(xyz中的一个)的栅格数；
         pad_count       number of cells beyond box

       Outputs:
         gridCenters  gridCounts(0)*gridCounts(1)*gridCounts(2) by 3 list of cell center positions
         gridCounts  三个维度上栅格的数目；
    */
  template <
    typename DerivedV,
    typename DerivedGV,
    typename Derivedside>
  IGL_INLINE void voxel_grid(
    const Eigen::MatrixBase<DerivedV> & vers, 
    const typename DerivedV::Scalar offset,
    const int s,
    const int pad_count,
    Eigen::PlainObjectBase<DerivedGV> & gridCenters,
    Eigen::PlainObjectBase<Derivedside> & gridCounts);
}
#ifndef IGL_STATIC_LIBRARY
#  include "voxel_grid.cpp"
#endif
#endif
