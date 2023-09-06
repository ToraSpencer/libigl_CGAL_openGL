#include <igl/cat.h>
#include <igl/edge_lengths.h>
#include <igl/parula.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/readMESH.h>
#include <igl/signed_distance.h>
#include <igl/slice_mask.h>
#include <igl/marching_tets.h>
#include <igl/upsample.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/writeOBJ.h>
#include <Eigen/Sparse>
#include <iostream>
#include "tutorial_shared_path.h"

// 符号距离场(signed distance field)

Eigen::MatrixXd vers;
Eigen::MatrixXi tets, tris;

igl::AABB<Eigen::MatrixXd,3> tree;          
igl::FastWindingNumberBVH fwn_bvh;

Eigen::MatrixXd normals,VN,EN;
Eigen::MatrixXi E;
Eigen::VectorXi edgeUeInfo;
double max_distance = 1;

double slice_z = 0.5;
bool overlay = false;

bool useFastWindingNumber = false;


// 每次键盘事件触发后更新viewer中的数据
void update_visualization(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;

  //    横截面网格：
  MatrixXd versSection;
  MatrixXi trisSection;

  // 1. 构造探测平面plane;
  Eigen::Vector4d plane(0,0,1,-((1-slice_z)*vers.col(2).minCoeff()+slice_z*vers.col(2).maxCoeff()));

  // 2. 选取探测平面切割的三角片；Extract triangle mesh slice through volume mesh and subdivide nasty triangles
  {
    VectorXi J;
    SparseMatrix<double> bary;

    {
      // Value of plane's implicit function at all vertices
      const VectorXd IV = (vers.col(0)*plane(0) +  vers.col(1)*plane(1) + vers.col(2)*plane(2)).array() + plane(3);
      igl::marching_tets(vers, tets, IV, versSection, trisSection, J, bary);
      igl::writeOBJ("vis.obj", versSection, trisSection);
    }

    while(true)
    {
      MatrixXd l;
      igl::edge_lengths(versSection, trisSection, l);
      l /= (versSection.colwise().maxCoeff() - versSection.colwise().minCoeff()).norm();
      const double max_l = 0.03;
      if(l.maxCoeff()<max_l)
        break;

      Array<bool,Dynamic,1> bad = l.array().rowwise().maxCoeff() > max_l;
      MatrixXi F_vis_bad, F_vis_good;
      igl::slice_mask(trisSection,bad,1,F_vis_bad);
      igl::slice_mask(trisSection,(bad!=true).eval(),1,F_vis_good);
      igl::upsample(versSection,F_vis_bad);
      trisSection = igl::cat(1,F_vis_bad,F_vis_good);
    }
  }

  // 3. 计算符号距离场；
  VectorXd S_vis; 
  if (!useFastWindingNumber)    // 如果输入网格是water-tight的，则可以不使用缠绕数；使用缠绕数的方法鲁棒性强但是速度慢；
  {
    VectorXi I;
    MatrixXd N,C;
    signed_distance_pseudonormal(versSection,vers,tris,tree, normals, VN,EN,edgeUeInfo,S_vis,I,C,N);     // bunny网格是water-tight的；
  } 
  else 
    signed_distance_fast_winding_number(versSection, vers, tris, tree, fwn_bvh, S_vis);


  //        lambda —— 直接融合网格
  const auto & append_mesh = [&trisSection,&versSection](const Eigen::MatrixXd & vers, const Eigen::MatrixXi & tris,const RowVector3d & color)
  {
    trisSection.conservativeResize(trisSection.rows()+tris.rows(),3);
    trisSection.bottomRows(tris.rows()) = tris.array()+versSection.rows();
    versSection.conservativeResize(versSection.rows()+vers.rows(),3);
    versSection.bottomRows(vers.rows()) = vers;
  };

  if(overlay)
    append_mesh(vers, tris, RowVector3d(0.8,0.8,0.8));

  viewer.data().clear();
  viewer.data().set_mesh(versSection, trisSection);
  viewer.data().set_data(S_vis);
  viewer.core().lighting_factor = overlay;
}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
{
  switch(key)
  {
    default:
      return false;
    case ' ':
      overlay ^= true;
      break;
    case '.':
      slice_z = std::min(slice_z+0.01,0.99);
      break;
    case ',':
      slice_z = std::max(slice_z-0.01,0.01);
      break;
    case '1':
      useFastWindingNumber = true;
      break;
    case '2':
      useFastWindingNumber = false;
      break;
  }
  update_visualization(viewer);
  return true;
}


int originalMain(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  cout<<"Usage:"<<endl;
  cout<<"[space]  toggle showing surface."<<endl;
  cout<<"'.'/','  push back/pull forward slicing plane."<<endl;
  cout<< "1/2 toggle between fast winding number (1) and pseudonormal (2) signing. \n";
  cout<<endl;

  // Load mesh: (vers,tets) tet-mesh of convex hull, tris contains original surface triangles
  igl::readMESH(TUTORIAL_SHARED_PATH "/bunny.mesh",vers, tets, tris);    

  // Encapsulated call to point_mesh_squared_distance to determine bounds
  {
    VectorXd sqrD;
    VectorXi I;
    MatrixXd closestVers;
    igl::point_mesh_squared_distance(vers, vers, tris, sqrD, I, closestVers);
    max_distance = sqrt(sqrD.maxCoeff());
  }

  // Fast winding and Pseudo normal depend on differnt AABB trees... We initialize both here.

  // Pseudonormal setup...
 
  // Precompute signed distance AABB tree
  tree.init(vers, tris);

  // Precompute vertex,edge and face normals
  igl::per_face_normals(vers, tris, normals);
  igl::per_vertex_normals(vers, tris, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, normals, VN);
  igl::per_edge_normals(vers,tris,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,normals,EN,E,edgeUeInfo);

  // fast winding number setup (just init fwn bvh)
  igl::fast_winding_number(vers, tris, 2, fwn_bvh);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  update_visualization(viewer);
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.launch();

  return 0;
}
 

int main(int argc, char** argv) 
{
    // return originalMain(argc, argv);

    test0();

    return 0;
}
