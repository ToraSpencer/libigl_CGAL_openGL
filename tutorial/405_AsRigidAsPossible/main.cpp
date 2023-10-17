#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>

#include <igl/arap.h>                               // ARAP变形
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

#include "tutorial_shared_path.h"

// ARAP deformation(As Rigid As Possible Deformation)——尽可能刚性的变形；

using RotationList = std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >; 

const Eigen::RowVector3d sea_green(70./255., 252./255., 167./255.); 
Eigen::MatrixXd vers0, vers; 
Eigen::MatrixXi tris0; 
Eigen::VectorXi S, b; 
Eigen::RowVector3d mid; 
double anim_t = 0.0; 
double anim_t_dir = 0.03; 
igl::ARAPData arap_data; 


bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen; 
  using namespace std; 

    MatrixXd bc(b.size(), vers0.cols()); 
    for(int i = 0; i < b.size(); i++)
    {
      bc.row(i) = vers0.row(b(i)); 
      switch(S(b(i)))
      {
        case 0:
        {
          const double r = mid(0)*0.25; 
          bc(i, 0) += r*sin(0.5*anim_t*2.*igl::PI); 
          bc(i, 1) -= r+r*cos(igl::PI+0.5*anim_t*2.*igl::PI); 
          break; 
        }

        case 1:
        {
          const double r = mid(1)*0.15; 
          bc(i, 1) += r+r*cos(igl::PI+0.15*anim_t*2.*igl::PI); 
          bc(i, 2) -= r*sin(0.15*anim_t*2.*igl::PI); 
          break; 
        }

        case 2:
        {
          const double r = mid(1)*0.15; 
          bc(i, 2) += r+r*cos(igl::PI+0.35*anim_t*2.*igl::PI); 
          bc(i, 0) += r*sin(0.35*anim_t*2.*igl::PI); 
          break; 
        }

        default:
          break; 
      }
    }

    igl::arap_solve(bc, arap_data, vers); 
    viewer.data().set_vertices(vers); 
    viewer.data().compute_normals(); 
    if(viewer.core().is_animating) 
        anim_t += anim_t_dir;  

  return false; 
}


bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating; 
      return true; 
  }
  return false; 
}


int main(int argc, char *argv[])      
{
  using namespace Eigen; 
  using namespace std; 


  igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off", vers0, tris0); 

  vers = vers0; 
  igl::readDMAT(TUTORIAL_SHARED_PATH "/decimated-knight-selection.dmat", S); 

  // vertices in selection
  igl::colon<int>(0, vers0.rows()-1, b); 

  auto tmp1 = stable_partition(b.data(), b.data() + b.size(), \
      [](int i)->bool
      {
          return S(i) >= 0;
      });
  auto tmp2 = b.data();
  auto tmp = stable_partition(b.data(), b.data() + b.size(), \
      [](int i)->bool
      {
          return S(i) >= 0;
      }) - b.data();

  b.conservativeResize(tmp); 
  auto tmp3 = tmp1 + 1;
  
  // Centroid
  mid = 0.5*(vers0.colwise().maxCoeff() + vers0.colwise().minCoeff()); 

  // Precomputation
  arap_data.max_iter = 100; 
  igl::arap_precomputation(vers0, tris0, vers0.cols(), b, arap_data); 

  // Set color based on selection
  MatrixXd C(tris0.rows(), 3); 
  RowVector3d purple(80.0/255.0, 64.0/255.0, 255.0/255.0); 
  RowVector3d gold(255.0/255.0, 228.0/255.0, 58.0/255.0); 
  for(int f = 0; f<tris0.rows(); f++)
  {
    if( S(tris0(f, 0))>=0 && S(tris0(f, 1))>=0 && S(tris0(f, 2))>=0) 
        C.row(f) = purple; 
    else
        C.row(f) = gold; 
  }


  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer; 
  viewer.data().set_mesh(vers, tris0); 
  viewer.data().set_colors(C); 
  viewer.callback_pre_draw = &pre_draw; 
  viewer.callback_key_down = &key_down; 
  viewer.core().is_animating = false; 
  viewer.core().animation_max_fps = 30.; 

  std::cout << "Press [space] to toggle animation" << endl; 

  viewer.launch(); 
}
