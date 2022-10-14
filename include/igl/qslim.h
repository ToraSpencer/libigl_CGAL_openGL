#ifndef IGL_QSLIM_H
#define IGL_QSLIM_H
#include "igl_inline.h"
#include <Eigen/Core>


namespace igl
{
    // qslim()¡ª¡ª
    /*
       Decimate (simplify) a triangle mesh in nD according to the paper "Simplifying Surfaces 
                    with Color and Texture using Quadric Error Metrics" by [Garland and Heckbert, 1987] (technically a followup to qslim). 

       The mesh can have open boundaries but should be edge-manifold.
  
       Inputs:
         vers               #vers by dim list of vertex positions. Assumes that vertices with
                                       infinite coordinates are "points at infinity" being used to close up
                                       boundary edges with faces. This allows special subspace quadrice for
                                       boundary edges: There should never be more than one "point at
                                       infinity" in a single triangle.
         tris               #tris by 3 list of triangle indices into vers
         max_m        desired number of output faces

       Outputs:
         versOut              #versOut by dim list of output vertex posistions (can be same ref as vers)
         trisOut               #trisOut by 3 list of output face indices into versOut (can be same ref as tris)
         newOldTrisInfo              #trisOut list of indices into tris of birth face
         newOldVersInfo               #versOut list of indices into vers of birth vertices
  */

  IGL_INLINE bool qslim(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const size_t max_m,
    Eigen::MatrixXd & versOut,
    Eigen::MatrixXi & trisOut,
    Eigen::VectorXi & newOldTrisInfo,
    Eigen::VectorXi & newOldVersInfo);
}


#ifndef IGL_STATIC_LIBRARY
#  include "qslim.cpp"
#endif
#endif
