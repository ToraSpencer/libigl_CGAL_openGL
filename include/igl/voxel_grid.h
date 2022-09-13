 
#ifndef IGL_VOXEL_GRID_H
#define IGL_VOXEL_GRID_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace igl
{

    // ������ά�����Χ��box��դ��
  /*
    Construct the cell center positions of a regular voxel grid(lattice) made  of perfectly square voxels.
   
       Inputs:
         box        bounding box to enclose by grid
         s          number of cell centers on largest gridCounts (including 2 * pad_count)
                    ��������Ǹ�ά��(xyz�е�һ��)��դ������
         pad_count       number of cells beyond box

       Outputs:
         gridCenters  gridCounts(0)*gridCounts(1)*gridCounts(2) by 3 list of cell center positions
         gridCounts  ����ά����դ�����Ŀ��
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


  // �����������vers�İ�Χ�е�դ��
    /*
 
       Inputs:
         vers        �������
         offset     ��Χ�б�Ե��������С���룻
         s          number of cell centers on largest gridCounts (including 2 * pad_count)
                    ��������Ǹ�ά��(xyz�е�һ��)��դ������
         pad_count       number of cells beyond box

       Outputs:
         gridCenters  gridCounts(0)*gridCounts(1)*gridCounts(2) by 3 list of cell center positions
         gridCounts  ����ά����դ�����Ŀ��
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
