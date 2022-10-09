#include "oriented_facets.h"

template <typename DerivedF, typename DerivedE>
IGL_INLINE void igl::oriented_facets(
  const Eigen::MatrixBase<DerivedF> & tris,
  Eigen::PlainObjectBase<DerivedE> & E)
{
  E.resize(tris.rows()*tris.cols(),tris.cols()-1);
  typedef typename DerivedE::Scalar EScalar;
  switch(tris.cols())
  {
    case 4:
      E.block(0*tris.rows(),0,tris.rows(),1) = tris.col(1).template cast<EScalar>();
      E.block(0*tris.rows(),1,tris.rows(),1) = tris.col(3).template cast<EScalar>();
      E.block(0*tris.rows(),2,tris.rows(),1) = tris.col(2).template cast<EScalar>();

      E.block(1*tris.rows(),0,tris.rows(),1) = tris.col(0).template cast<EScalar>();
      E.block(1*tris.rows(),1,tris.rows(),1) = tris.col(2).template cast<EScalar>();
      E.block(1*tris.rows(),2,tris.rows(),1) = tris.col(3).template cast<EScalar>();

      E.block(2*tris.rows(),0,tris.rows(),1) = tris.col(0).template cast<EScalar>();
      E.block(2*tris.rows(),1,tris.rows(),1) = tris.col(3).template cast<EScalar>();
      E.block(2*tris.rows(),2,tris.rows(),1) = tris.col(1).template cast<EScalar>();

      E.block(3*tris.rows(),0,tris.rows(),1) = tris.col(0).template cast<EScalar>();
      E.block(3*tris.rows(),1,tris.rows(),1) = tris.col(1).template cast<EScalar>();
      E.block(3*tris.rows(),2,tris.rows(),1) = tris.col(2).template cast<EScalar>();
      return;
    case 3:
      E.block(0*tris.rows(),0,tris.rows(),1) = tris.col(1).template cast<EScalar>();
      E.block(0*tris.rows(),1,tris.rows(),1) = tris.col(2).template cast<EScalar>();
      E.block(1*tris.rows(),0,tris.rows(),1) = tris.col(2).template cast<EScalar>();
      E.block(1*tris.rows(),1,tris.rows(),1) = tris.col(0).template cast<EScalar>();
      E.block(2*tris.rows(),0,tris.rows(),1) = tris.col(0).template cast<EScalar>();
      E.block(2*tris.rows(),1,tris.rows(),1) = tris.col(1).template cast<EScalar>();
      return;
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::oriented_facets<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1>, Eigen::Matrix<unsigned int, -1, 2, 0, -1, 2> >(Eigen::MatrixBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, 2, 0, -1, 2> >&);
// generated by autoexplicit.sh
template void igl::oriented_facets<Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::oriented_facets<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::oriented_facets<Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&);
template void igl::oriented_facets<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 2, 0, -1, 2> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> >&);
template void igl::oriented_facets<Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::oriented_facets<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 2, 0, -1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 2, 0, -1, 2> >&);
#endif
