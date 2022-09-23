#include "voxel_grid.h"
#include "grid.h"

template <
  typename Scalar, 
  typename DerivedGV, 
  typename Derivedside>
IGL_INLINE void igl::voxel_grid(
  const Eigen::AlignedBox<Scalar, 3> & box, 
  const int largestCount, 
  const int pad_count, 
  Eigen::PlainObjectBase<DerivedGV> & gridCenters, 
  Eigen::PlainObjectBase<Derivedside> & gridCounts)
{
  using namespace Eigen;
  using namespace std;
  typename DerivedGV::Index maxCompIdx = -1;            // 包围盒box的对角线向量中最大的分量的索引；0,1,2分别对应xyz分量；
  gridCounts.resize(1, 3);
  box.diagonal().maxCoeff(&maxCompIdx);                     
  const Scalar maxComp = box.diagonal()(maxCompIdx);          // 包围盒box的对角线向量中最大的分量；
 
  assert(largestCount > (pad_count*2+1) && "largestCount should be > 2*pad_count+1");

  // 计算xyz三个维度上栅格的个数gridCounts
  const Scalar largestCount0 = largestCount - 2*pad_count;
  gridCounts(maxCompIdx) = largestCount0;
  for(int i = 0; i < 3 ; i++)
  {
    if(i!=maxCompIdx)
        gridCounts(i) = std::ceil(largestCount0 * (box.diagonal()(i)) / maxComp);
        // gridCounts(i) = std::ceil(largestCount0 * (box.max()(i)-box.min()(i))/maxComp);
  }
  gridCounts.array() += 2 * pad_count;

  // 计算gridCenters;
  grid(gridCounts, gridCenters);

  /*
       A *    p/largestCount  + B = min
       A * (1-p/largestCount) + B = max
       B = min - A * p/largestCount
       A * (1-p/largestCount) + min - A * p/largestCount = max
       A * (1-p/largestCount) - A * p/largestCount = max-min
       A * (1-2p/largestCount) = max-min
       A  = (max-min)/(1-2p/largestCount)
  */

  const Array<Scalar, 3, 1> ps= (Scalar)(pad_count)/(gridCounts.transpose().template cast<Scalar>().array()-1.);
  const Array<Scalar, 3, 1> A = box.diagonal().array()/(1.0-2.*ps);

  /*
      // This would result in an "anamorphic",  but perfectly fit grid:
      const Array<Scalar, 3, 1> B = box.min().array() - A.array()*ps;
      gridCenters.array().rowwise() *= A.transpose();
      gridCenters.array().rowwise() += B.transpose();
       Instead scale by largest factor and move to match center
  */
  typename Array<Scalar, 3, 1>::Index ai = -1;
  Scalar a = A.maxCoeff(&ai);
  const Array<Scalar, 1, 3> ratio =  a*(gridCounts.template cast<Scalar>().array()-1.0)/(Scalar)(gridCounts(ai)-1.0);
  gridCenters.array().rowwise() *= ratio;
  const Eigen::Matrix<Scalar, 1, 3> offset = (box.center().transpose()-gridCenters.colwise().mean()).eval();
  gridCenters.rowwise() += offset;
}



template <
  typename DerivedV, 
  typename DerivedGV, 
  typename Derivedside>
IGL_INLINE void igl::voxel_grid(
  const Eigen::MatrixBase<DerivedV> & V,  
  const typename DerivedV::Scalar offset, 
  const int largestCount, 
  const int pad_count, 
  Eigen::PlainObjectBase<DerivedGV> & gridCenters, 
  Eigen::PlainObjectBase<Derivedside> & gridCounts)
{
  typedef typename DerivedV::Scalar Scalar;
  Eigen::AlignedBox<Scalar, 3> box;
  typedef Eigen::Matrix<Scalar, 1, 3> RowVector3S;
  assert(V.cols() == 3 && "V must contain positions in 3D");
  RowVector3S min_ext = V.colwise().minCoeff().array() - offset;
  RowVector3S max_ext = V.colwise().maxCoeff().array() + offset;
  box.extend(min_ext.transpose());
  box.extend(max_ext.transpose());
  return igl::voxel_grid(box, largestCount, 1, gridCenters, gridCounts);
}


// 模板特化：
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::voxel_grid<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>::Scalar,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);
// generated by autoexplicit.sh
template void igl::voxel_grid<Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>::Scalar,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);
// generated by autoexplicit.sh
template void igl::voxel_grid<Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> > const&,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>::Scalar,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);
// generated by autoexplicit.sh
template void igl::voxel_grid<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>::Scalar,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);
// generated by autoexplicit.sh
template void igl::voxel_grid<float,  Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  3,  1,  0,  3,  1> >(Eigen::AlignedBox<float,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  3,  1,  0,  3,  1> >&);
// generated by autoexplicit.sh
template void igl::voxel_grid<float,  Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::AlignedBox<float,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);
template void igl::voxel_grid<float,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  3,  1,  0,  3,  1> >(Eigen::AlignedBox<float,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  3,  1,  0,  3,  1> >&);
template void igl::voxel_grid<float,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  3,  1,  0,  3,  1> >(Eigen::AlignedBox<float,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  3,  1,  0,  3,  1> >&);
template void igl::voxel_grid<double,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  3,  1,  0,  3,  1> >(Eigen::AlignedBox<double,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  3,  1,  0,  3,  1> >&);
template void igl::voxel_grid<double,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  1,  3,  1,  1,  3> >(Eigen::AlignedBox<double,  3> const&,  int,  int,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  1,  3,  1,  1,  3> >&);



#endif
