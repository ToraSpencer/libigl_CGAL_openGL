#include <igl/copyleft/marching_cubes.h>
#include <igl/unique_simplices.h>
#include <igl/dual_contouring.h>
#include <igl/get_seconds.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <igl/per_corner_normals.h>
#include <igl/parallel_for.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/polygon_corners.h>
#include <igl/slice.h>
#include <igl/sparse_voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>


int main(int argc,  char * argv[])
{
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds(); 
    double diff = igl::get_seconds()-t_start; 
    t_start += diff; 
    return diff; 
  }; 

  // 1. Create an interesting shape with sharp features using SDF CSG with spheres.
  const auto& sphere = [](const Eigen::RowVector3d & c,  const double r,   const Eigen::RowVector3d & x)->double
  {
    return (x-c).norm() - r; 
  }; 

  const std::function<double(const Eigen::RowVector3d & x)> tris0 =  [&](const Eigen::RowVector3d & x)->double
  {
    return 
      std::min( std::min(std::max(
          sphere(Eigen::RowVector3d(-0.2, 0, -0.2), 0.5, x), 
          -sphere(Eigen::RowVector3d(+0.2, 0, 0.2), 0.5, x)), 
           sphere(Eigen::RowVector3d(-0.15, 0, -0.15), 0.3, x)
         ), 
        std::max( std::max( sphere(Eigen::RowVector3d(-0.2, -0.5, -0.2), 0.6, x), x(1)+0.45), -0.6-x(1)) ); 
  }; 
  Eigen::RowVector3d p0(-0.2, 0.5, -0.2); 
  assert(abs(tris0(p0)) < 1e-10 && "p0 should be on zero level-set"); 

  // 2. Simple finite difference gradients
  const auto & fd = [](const std::function<double(const Eigen::RowVector3d&)>& tris,  const Eigen::RowVector3d & x)
  {
    const double eps = 1e-10; 
    Eigen::RowVector3d g; 
    for(int c = 0; c<3; c++)
    {
      const Eigen::RowVector3d xp = x+eps*Eigen::RowVector3d(c==0, c==1, c==2); 
      const double fp = tris(xp); 
      const Eigen::RowVector3d xn = x-eps*Eigen::RowVector3d(c==0, c==1, c==2); 
      const double fn = tris(xn); 
      g(c) = (fp-fn)/(2*eps); 
    }
    return g; 
  }; 

  const auto & f_grad = [&fd, &tris0](const Eigen::RowVector3d & x)
  {
    return fd(tris0, x).normalized(); 
  }; 

  Eigen::MatrixXd vers; 
  Eigen::MatrixXi Q, tris; 
  Eigen::MatrixXd versResult_SDF, mcN; 
  Eigen::MatrixXi trisResult_SDF; 

  // Grid parameters
  const Eigen::RowVector3d min_corner(-2, -2, -2); 
  const Eigen::RowVector3d max_corner(+2, +2, +2); 
  const int s = 256; 
  int gridCountsX = s+1; 
  int gridCountsY = s+1; 
  int gridCountsZ = s+1; 
  const Eigen::RowVector3d step =  (max_corner-min_corner).array()/(Eigen::RowVector3d(gridCountsX, gridCountsY, gridCountsZ).array()-1); 

  // Sparse grid below assumes regular grid
  assert((step(0) == step(1))&&(step(0) == step(2))); 

  // Dual contouring parameters
  bool constrained = false; 
  bool triangles = false; 
  bool root_finding = true; 
  for(int pass = 0; pass < 2; pass++)
  {
    const bool sparse = pass == 1; 
    printf("Using %s grid..\n", sparse?"sparse":"dense"); 
    if(sparse)
    {
      // igl::sparse_voxel_grid assumes (0, 0, 0) lies on the grid. But dense igl::grid
      // below won't necessarily do that depending on gridCountsX, gridCountsY, gridCountsZ.
      tictoc(); 
      Eigen::MatrixXd gridCenters; 
      Eigen::VectorXd SDF; 
      Eigen::Matrix<int, Eigen::Dynamic, 8> GI; 
      igl::sparse_voxel_grid(p0, tris0, step(0), 16.*pow(step(0), -2.), SDF, gridCenters, GI); 
      const auto t_Gf = tictoc(); 
      printf("  %5f secs to populate sparse grid of %ld cells\n", t_Gf+tictoc(), GI.rows()); 

      // Dual contouring requires list of sparse edges (not cells), extract _all_ edges from sparse_voxel_grid (conservative)
      Eigen::Matrix<int, Eigen::Dynamic, 2> GI2; 
      {
        Eigen::Matrix<int, Eigen::Dynamic, 2> all_GI2(GI.rows()*12, 2); 
        all_GI2 << 
          // front
          GI.col(0), GI.col(1), 
          GI.col(1), GI.col(2), 
          GI.col(2), GI.col(3), 
          GI.col(3), GI.col(0), 
          // back
          GI.col(4+0), GI.col(4+1), 
          GI.col(4+1), GI.col(4+2), 
          GI.col(4+2), GI.col(4+3), 
          GI.col(4+3), GI.col(4+0), 
          // sides
          GI.col(0), GI.col(4+0),  
          GI.col(1), GI.col(4+1),  
          GI.col(2), GI.col(4+2), 
          GI.col(3), GI.col(4+3); 
        Eigen::VectorXi _1, _2; 
        igl::unique_simplices(all_GI2, GI2, _1, _2); 
      }

      // dual contouring
      tictoc(); 
      Eigen::RowVector3d step = (max_corner-min_corner).array()/(Eigen::RowVector3d(gridCountsX, gridCountsY, gridCountsZ).array()-1); 
      igl::dual_contouring(tris0, f_grad, step, SDF, gridCenters, GI2, constrained, triangles, root_finding, vers, Q); 
      printf("  %5f secs dual contouring\n", t_Gf+tictoc());

      // marching cubes, Could use igl::marching_cubes once https://github.com/libigl/libigl/pull/1687 is merged
      tictoc(); 
      igl::copyleft::marching_cubes(SDF, gridCenters, GI, versResult_SDF,  trisResult_SDF); 
      printf("  %5f secs marching cubes\n", t_Gf+tictoc()); 
    }
    else
    {
      tictoc(); 
      igl::dual_contouring(tris0, f_grad, min_corner, max_corner, gridCountsX, gridCountsY, gridCountsZ, constrained, triangles, root_finding, vers, Q); 
      printf("  %5f secs dual contouring\n", tictoc()); 

      // build and sample grid
      tictoc(); 
      Eigen::MatrixXd gridCenters; 
      igl::grid(Eigen::RowVector3i(gridCountsX, gridCountsY, gridCountsZ), gridCenters); 
      Eigen::VectorXd SDF(gridCenters.rows()); 
      igl::parallel_for(gridCenters.rows(), [&](const int i)
      {
        gridCenters.row(i).array() *= (max_corner-min_corner).array(); 
        gridCenters.row(i) += min_corner; 
        SDF(i) = tris0(gridCenters.row(i)); 
      }, 1000ul); 
      const auto t_grid = tictoc(); 
      igl::marching_cubes(SDF, gridCenters, gridCountsX, gridCountsY, gridCountsZ, 0, versResult_SDF, trisResult_SDF); 
      const auto t_mc = tictoc(); 
      printf("  %5f secs (%5f + %5f) marching cubes\n", t_grid+t_mc, t_grid, t_mc); 
    }
  }

  // Crisp (as possible) rendering of resulting MC triangle mesh
  igl::per_corner_normals(versResult_SDF, trisResult_SDF, 20, mcN); 

  // Crisp rendering of resulting DC quad mesh with edges
  Eigen::MatrixXi edges; 
  Eigen::MatrixXd VV, N, NN; 
  Eigen::VectorXi J; 
  Eigen::MatrixXi FF; 
  if(triangles)
  {
    VV = vers; 
    FF = Q; 
    edges.resize(Q.rows()*3, 2); 
    edges<<
      Q.col(0),  Q.col(1),  
      Q.col(1),  Q.col(2),  
      Q.col(2),  Q.col(0); 
  }
  else
  {
    Eigen::VectorXi I, C; 
    igl::polygon_corners(Q, I, C); 
    edges.resize(Q.rows()*4, 2); 
    edges<<
      Q.col(0),  Q.col(1),  
      Q.col(1),  Q.col(2),  
      Q.col(2),  Q.col(3),  
      Q.col(3),  Q.col(0); 
    igl::per_face_normals(vers, I, C, N, VV, FF, J); 
    igl::slice(N, J, 1, NN); 
    igl::per_corner_normals(vers, I, C, 20, N, VV, FF, J, NN); 
  }


  // visualize:
  igl::opengl::glfw::Viewer viewer; 
  bool show_edges = true; 
  bool use_dc = true; 

  const auto update = [&]()
  {
    const bool was_face_based = viewer.data().face_based ; 
    viewer.data().clear(); 
    if(use_dc)
    {
      viewer.data().set_mesh(VV, FF); 
      viewer.data().show_lines = false; 
      viewer.data().set_normals(NN); 
      if(show_edges)
      {
        viewer.data().clear_edges(); 
        viewer.data().set_edges(vers, edges, Eigen::RowVector3d(0, 0, 0)); 
      }
    }
    else
    {
      viewer.data().set_mesh(versResult_SDF, trisResult_SDF); 
      viewer.data().set_normals(mcN); 
      viewer.data().show_lines = show_edges; 
    }
    viewer.data().face_based = was_face_based; 
  }; 
  update(); 
  viewer.data().face_based = true; 
  viewer.callback_key_pressed = [&](decltype(viewer) &, unsigned int key,  int mod)
  {
    switch(key)
    {
      case ' ': use_dc=!use_dc;  update(); return true; 
      case 'L': case 'l': show_edges=!show_edges;  update(); return true; 
    }
    return false; 
  }; 

  std::cout<<R"([space]  Toggle between dual contouring and marching cubes)"; 
  viewer.launch(); 
}

