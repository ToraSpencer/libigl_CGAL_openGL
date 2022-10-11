#include "collapse_small_triangles.h"

#include "bounding_box_diagonal.h"
#include "doublearea.h"
#include "edge_lengths.h"
#include "colon.h"
#include "faces_first.h"

#include <limits>
#include <iostream>

void igl::collapse_small_triangles(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const double eps, 
  Eigen::MatrixXi & trisOut)
{
  using namespace Eigen; 
  using namespace std; 

  // Compute bounding box diagonal length
  double boxArrow = bounding_box_diagonal(vers);            // 包围盒的对角线长度；
  MatrixXd triEdgeLen; 
  edge_lengths(vers,  tris,  triEdgeLen); 
  VectorXd dbArea;                                                      // 三角片面积的两倍；
  doublearea(triEdgeLen, 0., dbArea); 

  // 判定三角片是否要折叠的面积阈值；
  const double min_dblarea = 2.0 * eps * boxArrow * boxArrow;       

  Eigen::VectorXi FIM = colon<int>(0, vers.rows()-1); 
  int num_edge_collapses = 0; 

  // Loop over triangles
  for(int triIdx = 0; triIdx < tris.rows(); triIdx++)
  {
    if(dbArea(triIdx) < min_dblarea)
    {
      double minl = 0; 
      int minli = -1; 

      // Find shortest edge
      for(int e = 0; e < 3; e++)
      {
        if(minli==-1 || triEdgeLen(triIdx, e)<minl)
        {
          minli = e; 
          minl = triEdgeLen(triIdx, e); 
        }
      }

      double maxl = 0; 
      int maxli = -1; 

      // Find longest edge
      for(int e = 0; e<3; e++)
      {
        if(maxli==-1 || triEdgeLen(triIdx, e)>maxl)
        {
          maxli = e; 
          maxl = triEdgeLen(triIdx, e); 
        }
      }

      // Be sure that min and max aren't the same
      maxli = (minli==maxli?(minli+1)%3:maxli); 

      // Collapse min edge maintaining max edge: i-->j;  Q: Why this direction?
      int i = maxli; 
      int j = ((minli+1)%3 == maxli ? (minli+2)%3: (minli+1)%3); 
      assert(i != minli); 
      assert(j != minli); 
      assert(i != j); 
      FIM(tris(triIdx, i)) = FIM(tris(triIdx, j)); 
      num_edge_collapses++; 
    }
  }

  // Reindex faces
  MatrixXi tris0 = tris; 

  // Loop over triangles
  for(int triIdx = 0; triIdx<tris0.rows(); triIdx++)
        for(int i = 0; i<tris0.cols(); i++)
            tris0(triIdx, i) = FIM(tris0(triIdx, i)); 

  trisOut.resizeLike(tris0); 
  int num_face_collapses=0; 

  // Only keep uncollapsed faces
  {
    int index = 0; 

    // Loop over triangles
    for(int triIdx = 0; triIdx<tris0.rows(); triIdx++)
    {
      bool collapsed = false; 

      // Check if any indices are the same
      for(int i = 0; i<tris0.cols(); i++)
      {
        for(int j = i+1; j < tris0.cols(); j++)
        {
          if(tris0(triIdx, i)==tris0(triIdx, j))
          {
            collapsed = true; 
            num_face_collapses++; 
            break; 
          }
        }
      }
      if(!collapsed)
            trisOut.row(index++) = tris0.row(triIdx); 
    }

    // Use conservative resize
    trisOut.conservativeResize(index, trisOut.cols()); 
  }

  if(num_edge_collapses == 0)
  {
    // There must have been a "collapsed edge" in the input
    assert(num_face_collapses==0); 

    // Base case
    return; 
  }

  MatrixXi recFF = trisOut; 

  return collapse_small_triangles(vers, recFF, eps, trisOut); 
}
