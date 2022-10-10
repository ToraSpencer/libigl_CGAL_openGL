#include "barycenter.h"


// 计算网格每个三角片的重心坐标
template <  typename DerivedV,   typename DerivedF,   typename DerivedBC>
IGL_INLINE void igl::barycenter( const Eigen::MatrixBase<DerivedV> & vers,    \
            const Eigen::MatrixBase<DerivedF> & tris, \
             Eigen::PlainObjectBase<DerivedBC> & baryCenters)
{
    baryCenters.setZero(tris.rows(), vers.cols()); 

  for(int i = 0; i<tris.rows(); i++)        // Loop over faces
  {
    // loop around face
    for(int j = 0; j<tris.cols(); j++)
        baryCenters.row(i) += vers.row(tris(i, j));             // Accumulate
 
    baryCenters.row(i) /= double(tris.cols());       // average
  }
}



#ifdef IGL_STATIC_LIBRARY

template void igl::barycenter<Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&); 
template void igl::barycenter<Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<float,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >&); 
template void igl::barycenter<Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  4,  0,  -1,  4> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  4,  0,  -1,  4> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  4,  0,  -1,  4>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  4,  0,  -1,  4> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  4,  0,  -1,  4> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  4,  0,  -1,  4> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  0,  -1,  3>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  3,  0,  -1,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  2,  3,  0,  2,  3>,  Eigen::Matrix<double,  2,  3,  0,  2,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<double,  2,  3,  0,  2,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  2,  3,  0,  2,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  2,  3,  0,  2,  3>,  Eigen::Matrix<double,  2,  3,  0,  2,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  2,  3,  0,  2,  3> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  2,  3,  0,  2,  3> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  2,  0,  -1,  2> >&); 
template void igl::barycenter<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1>,  Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::MatrixBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  0,  -1,  3> >&); 
#endif
