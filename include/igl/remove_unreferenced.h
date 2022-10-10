#ifndef IGL_REMOVE_UNREFERENCED_H
#define IGL_REMOVE_UNREFERENCED_H
#include "igl_inline.h"
#include <Eigen/Core>


namespace igl 
{
    // remove_unreferenced()——去除网格中的孤立点；

    // 重载1.1.1
    /*
       Remove unreferenced vertices from vers, updating tris accordingly
  
       Input:
         vers                   #vers by dim list of mesh vertex positions
         tris                   #tris by ss list of simplices (Values of -1 are quitely skipped)

       Outputs:
         versOut              list of mesh vertex positions
         trisOut               #trisOut by ss list of simplices
         I                  #vers by 1 list of indices such that: trisOut = IM(tris) and NT = IM(T) and vers(find(IM<=size(versOut,1)),:) = versOut
         J                 #versOut by 1 list, such that versOut = vers(J, :)
    */     
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI>
  IGL_INLINE void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &vers,
    const Eigen::MatrixBase<DerivedF> &tris,
    Eigen::PlainObjectBase<DerivedNV> &versOut,
    Eigen::PlainObjectBase<DerivedNF> &trisOut,
    Eigen::PlainObjectBase<DerivedI> &I);


  // 重载1.1
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI,
    typename DerivedJ>
  IGL_INLINE void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &vers,
    const Eigen::MatrixBase<DerivedF> &tris,
    Eigen::PlainObjectBase<DerivedNV> &versOut,
    Eigen::PlainObjectBase<DerivedNF> &trisOut,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);


  // 重载1.
  /*
       Inputs:
         versCount                  number of vertices (possibly greater than tris.maxCoeff()+1)
         tris               #tris by ss list of simplices

       Outputs:
         IM                 #vers by 1 list of indices such that: trisOut = IM(tris) and NT = IM(T)
                                            and vers(find(IM<=size(versOut,1)),:) = versOut
         J                      #RV by 1 list, such that RV = vers(J,:)
  */   
  template <
    typename DerivedF,
    typename DerivedI,
    typename DerivedJ>
  IGL_INLINE void remove_unreferenced(
    const size_t versCount,
    const Eigen::MatrixBase<DerivedF> &tris,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J);

}

#ifndef IGL_STATIC_LIBRARY
#  include "remove_unreferenced.cpp"
#endif

#endif
