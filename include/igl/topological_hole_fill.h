#ifndef IGL_TOPOLOGICAL_HOLE_FILL_H
#define IGL_TOPOLOGICAL_HOLE_FILL_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
    // topological_hole_fill()����������

    /*
           Topological fill hole on a mesh, with one additional vertex each hole

           Index of new abstract vertices will be tris.maxCoeff() + (index of hole)
  
           Inputs:
             tris                  #tris by simplex-size list of element indices
             b                    #b boundary indices to preserve
             holes              ˫��std::vector����һ�㣺ÿ��Ԫ�ش���һ�������ڶ��㣺ÿ��Ԫ�ش����Ķ���������
                                    vector of hole loops to fill

           Outputs:
             trisOut            input tris stacked with filled triangles.
  */
  template <
  typename DerivedF,
  typename Derivedb,
  typename VectorIndex,
  typename DerivedF_filled>
IGL_INLINE void topological_hole_fill(
  const Eigen::MatrixBase<DerivedF> & tris,
  const Eigen::MatrixBase<Derivedb> & b,
  const std::vector<VectorIndex> & holes,
  Eigen::PlainObjectBase<DerivedF_filled> &trisOut);
}


#ifndef IGL_STATIC_LIBRARY
#  include "topological_hole_fill.cpp"
#endif

#endif