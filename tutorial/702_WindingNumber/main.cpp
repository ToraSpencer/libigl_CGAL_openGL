#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/parula.h>
#include <igl/readMESH.h>
#include <igl/slice.h>
#include <igl/marching_tets.h>
#include <igl/winding_number.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>

#include "tutorial_shared_path.h"

Eigen::MatrixXd vers, baryCenters; 
Eigen::VectorXd windingNums; 
Eigen::MatrixXi tets, tris, G; 

double slice_z = 0.5;
enum OverLayType
{
    OVERLAY_NONE = 0,
    OVERLAY_INPUT = 1,
    OVERLAY_OUTPUT = 2,
    NUM_OVERLAY = 3,
} overlay = OVERLAY_NONE;


void update_visualization(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen; 
  using namespace std; 
  Eigen::Vector4d plane(0, 0, 1, -((1-slice_z)*vers.col(2).minCoeff()+slice_z*vers.col(2).maxCoeff())); 
  MatrixXd V_vis; 
  MatrixXi F_vis; 
  VectorXi J; 

  {
    SparseMatrix<double> bary; 
    // Value of plane's implicit function at all vertices
    const VectorXd IV = 
      (vers.col(0)*plane(0) + 
        vers.col(1)*plane(1) + 
        vers.col(2)*plane(2)).array()
      + plane(3); 
    igl::marching_tets(vers, tets, IV, V_vis, F_vis, J, bary); 
  }

  VectorXd W_vis; 
  igl::slice(windingNums, J, W_vis); 
  MatrixXd C_vis; 
  // color without normalizing
  igl::parula(W_vis, false, C_vis); 

  // lambda——
  const auto & append_mesh = [&C_vis, &F_vis, &V_vis](
        const Eigen::MatrixXd & vers, 
        const Eigen::MatrixXi & tris, 
        const RowVector3d & color)
  {
    F_vis.conservativeResize(F_vis.rows()+tris.rows(), 3); 
    F_vis.bottomRows(tris.rows()) = tris.array()+V_vis.rows(); 
    V_vis.conservativeResize(V_vis.rows()+vers.rows(), 3); 
    V_vis.bottomRows(vers.rows()) = vers; 
    C_vis.conservativeResize(C_vis.rows()+tris.rows(), 3); 
    C_vis.bottomRows(tris.rows()).rowwise() = color; 
  }; 


  switch(overlay)
  {
    case OVERLAY_INPUT:
      append_mesh(vers, tris, RowVector3d(1., 0.894, 0.227)); 
      break; 
    case OVERLAY_OUTPUT:
      append_mesh(vers, G, RowVector3d(0.8, 0.8, 0.8)); 
      break; 
    default:
      break; 
  }

  viewer.data().clear(); 
  viewer.data().set_mesh(V_vis, F_vis); 
  viewer.data().set_colors(C_vis); 
  viewer.data().set_face_based(true); 
}


bool key_down(igl::opengl::glfw::Viewer& viewer,  unsigned char key,  int mod)
{
  switch(key)
  {
    default:
      return false; 
    case ' ':
      overlay = (OverLayType)((1+(int)overlay)%NUM_OVERLAY); 
      break; 
    case '.':
      slice_z = std::min(slice_z+0.01, 0.99); 
      break; 
    case ', ':
      slice_z = std::max(slice_z-0.01, 0.01); 
      break; 
  }
  update_visualization(viewer); 
  return true; 
}


int main(int argc,  char *argv[])
{
  using namespace Eigen; 
  using namespace std; 

  cout<<"Usage:"<<endl; 
  cout<<"[space]  toggle showing input mesh,  output mesh or slice "<<endl; 
  cout<<"         through tet-mesh of convex hull."<<endl; 
  cout<<"'.'/', '  push back/pull forward slicing plane."<<endl; 
  cout<<endl; 

  // Load mesh: (vers, tets) tet-mesh of convex hull,  tris contains facets of input surface mesh _after_ self-intersection resolution
  igl::readMESH(TUTORIAL_SHARED_PATH "/big-sigcat.mesh", vers,  tets,  tris); 

  igl::barycenter(vers, tets, baryCenters);         // 计算每个体素的重心；

  // 计算每个体素重心的缠绕数；
  cout<<"Computing winding number over all "<< tets.rows()<< " tets..." << endl; 
  igl::winding_number(vers, tris, baryCenters, windingNums); 

  // Extract interior tets
  MatrixXi CT((windingNums.array()>0.5).count(), 4); 
  {
    size_t k = 0; 
    for(size_t t = 0; t<tets.rows(); t++)
    {
      if(windingNums(t)>0.5)
      {
        CT.row(k) = tets.row(t); 
        k++; 
      }
    }
  }

  // find bounary facets of interior tets
  igl::boundary_facets(CT, G); 

  // boundary_facets seems to be reversed...
  G = G.rowwise().reverse().eval(); 

  // normalize
  windingNums = (windingNums.array() - windingNums.minCoeff())/(windingNums.maxCoeff()-windingNums.minCoeff()); 

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer; 
  update_visualization(viewer); 
  viewer.callback_key_down = &key_down; 
  viewer.launch(); 
}
