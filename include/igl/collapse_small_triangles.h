#ifndef IGL_COLLAPSE_SMALL_TRIANGLES_H
#define IGL_COLLAPSE_SMALL_TRIANGLES_H
#include <Eigen/Dense>


namespace igl
{
    // collapse_small_triangles()¡ª¡ª

      /*
           Given a triangle mesh (vers,tris) compute a new mesh (VV, trisOut) which contains the original faces and vertices of (vers, tris) 
                except any small triangles have been removed via collapse.
  
           We are *not* following the rules in "Mesh Optimization" [Hoppe et al] Section 4.2. 
                But for our purposes we don't care about this criteria.
  
           Inputs:
                 vers           #vers by 3 list of vertex positions
                 tris            #tris by 3 list of triangle indices into vers
                 eps           epsilon for smallest allowed area treated as fraction of squared bounding box diagonal

           Outputs:
                 trisOut             #trisOut by 3 list of triangle indices into vers
  
      */


  void collapse_small_triangles(
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const double eps,
    Eigen::MatrixXi & trisOut);
}

#ifndef IGL_STATIC_LIBRARY
#  include "collapse_small_triangles.cpp"
#endif
  
#endif
