#include "connect_boundary_to_infinity.h"
#include "boundary_facets.h"

template <typename DerivedF,  typename DerivedFO>
IGL_INLINE void igl::connect_boundary_to_infinity(
  const Eigen::MatrixBase<DerivedF> & tris, 
  Eigen::PlainObjectBase<DerivedFO> & FO)
{
  return connect_boundary_to_infinity(tris, tris.maxCoeff(), FO);
}


template <typename DerivedF,  typename DerivedFO>
IGL_INLINE void igl::connect_boundary_to_infinity(
  const Eigen::MatrixBase<DerivedF> & tris, 
  const typename DerivedF::Scalar inf_index, 
  Eigen::PlainObjectBase<DerivedFO> & FO)
{
  // Determine boundary edges
  Eigen::Matrix<typename DerivedFO::Scalar, Eigen::Dynamic, Eigen::Dynamic> O;
  boundary_facets(tris, O);

  FO.resize(tris.rows()+O.rows(), tris.cols());
  typedef Eigen::Matrix<typename DerivedFO::Scalar, Eigen::Dynamic, 1> VectorXI;
  FO.topLeftCorner(tris.rows(), tris.cols()) = tris;
  FO.bottomLeftCorner(O.rows(), O.cols()) = O.rowwise().reverse();
  FO.bottomRightCorner(O.rows(), 1).setConstant(inf_index);
}


template <
  typename DerivedV,  
  typename DerivedF,  
  typename DerivedVO,  
  typename DerivedFO>
IGL_INLINE void igl::connect_boundary_to_infinity(
  const Eigen::MatrixBase<DerivedV> & vers, 
  const Eigen::MatrixBase<DerivedF> & tris, 
  Eigen::PlainObjectBase<DerivedVO> & VO, 
  Eigen::PlainObjectBase<DerivedFO> & FO)
{
  typename DerivedV::Index inf_index = vers.rows();
  connect_boundary_to_infinity(tris, inf_index, FO);
  VO.resize(vers.rows()+1, vers.cols());
  VO.topLeftCorner(vers.rows(), vers.cols()) = vers;
  auto inf = std::numeric_limits<typename DerivedVO::Scalar>::infinity();
  VO.row(vers.rows()).setConstant(inf);
}


#ifdef IGL_STATIC_LIBRARY
template void igl::connect_boundary_to_infinity<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&);
#endif
