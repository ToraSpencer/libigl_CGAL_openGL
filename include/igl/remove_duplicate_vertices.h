#ifndef IGL_REMOVE_DUPLICATE_VERTICES_H
#define IGL_REMOVE_DUPLICATE_VERTICES_H
#include "igl_inline.h"
#include <Eigen/Dense>


namespace igl
{

  // REMOVE_DUPLICATE_VERTICES Remove duplicate vertices upto a uniqueness tolerance (epsilon)
  /*
       Inputs:
         vers                #vers by dim list of vertex positions
         epsilon            uniqueness tolerance used coordinate-wise: 1e0 --> integer match, 1e-1 --> match up to first decimal, ... , 0 --> exact match.

       Outputs:
         versOut                    # versOut by dim new list of vertex positions
         selectedIdxes           # versOut by 1 list of indices so versOut = vers(selectedIdxes, :) 
         SVJ                           # vers by 1 list of indices so vers = versOut(SVJ, :)
  
       Example:
         % Mesh in (vers,tris)
         [versOut, selectedIdxes, SVJ] = remove_duplicate_vertices(vers,1e-7);

         % remap faces
         SF = SVJ(tris);
  */
  template <
    typename DerivedV, 
    typename DerivedSV, 
    typename DerivedSVI, 
    typename DerivedSVJ>
  IGL_INLINE void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& vers,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& versOut,
    Eigen::PlainObjectBase<DerivedSVI>& selectedIdxes,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ);


  // Wrapper that also remaps given faces (tris) --> (SF) so that SF index versOut
  template <
    typename DerivedV, 
    typename DerivedF,
    typename DerivedSV, 
    typename DerivedSVI, 
    typename DerivedSVJ,
    typename DerivedSF>
  IGL_INLINE void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& vers,
    const Eigen::MatrixBase<DerivedF>& tris,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& versOut,
    Eigen::PlainObjectBase<DerivedSVI>& selectedIdxes,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ,
    Eigen::PlainObjectBase<DerivedSF>& SF);
}

#ifndef IGL_STATIC_LIBRARY
#  include "remove_duplicate_vertices.cpp"
#endif

#endif
