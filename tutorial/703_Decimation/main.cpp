#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#include "tutorial_shared_path.h"


// ���񾫼򡪡����ڱߵľ���
int main(int argc,  char * argv[])
{
  using namespace std; 
  using namespace Eigen; 
  using namespace igl; 
  cout<<"Usage: ./703_Decimation_bin [filename.(off|obj|ply)]"<<endl; 
  cout<<"  [space]  toggle animation."<<endl; 
  cout<<"  'r'  reset."<<endl; 

  string filename(TUTORIAL_SHARED_PATH "/fertility.off"); 
  if(argc>=2)
        filename = argv[1]; 

  MatrixXd vers, versOld; 
  MatrixXi tris, trisOld; 
  igl::opengl::glfw::Viewer viewer; 
  read_triangle_mesh(filename, versOld, trisOld); 

  // Prepare array-based edge data structures and priority queue
  VectorXi EMAP; 
  MatrixXi edges, EF, EI; 
  igl::min_heap< std::tuple<double, int, int> > workQueue; 
  Eigen::VectorXi EQ; 

  // If an edge were collapsed,  we'd collapse it to these points:
  MatrixXd C; 
  int num_collapsed;            // ������ѭ���ļ�����

  // lambda�����������ݸ�λ
  const auto & reset = [&]()
  {
    tris = trisOld; 
    vers = versOld; 
    edge_flaps(tris, edges, EMAP, EF, EI); 
    C.resize(edges.rows(), vers.cols()); 
    VectorXd costs(edges.rows()); 

    // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
    workQueue = {}; 
    EQ = Eigen::VectorXi::Zero(edges.rows()); 

    {
      Eigen::VectorXd costs(edges.rows()); 
      igl::parallel_for(edges.rows(), \
          [&](const int e)
          {
            double cost = e; 
            RowVectorXd p(1, 3); 
            shortest_edge_and_midpoint(e, vers, tris, edges, EMAP, EF, EI, cost, p); 
            C.row(e) = p; 
            costs(e) = cost; 
          },\
          10000); 

      for(int e = 0; e<edges.rows(); e++)
        workQueue.emplace(costs(e), e, 0); 
    }

    num_collapsed = 0; 
    viewer.data().clear(); 
    viewer.data().set_mesh(vers, tris); 
    viewer.data().set_face_based(true); 
  }; 

  // lambda����ִ�����񾫼�׼����Ⱦ�����ݣ�
  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool
  {
    // ÿһ�ζ���ѭ���У�����10%�ıߣ�
    if(viewer.core().is_animating && !workQueue.empty())
    {
      bool FlagCollapsed = false;           // ����ѭ���б������Ƿ�ִ�гɹ���

      // ִ�б���������collapse edge
      const int max_iter = std::ceil(0.01*workQueue.size()); 
      for(int j = 0; j < max_iter; j++)
      {
        if(!collapse_edge(shortest_edge_and_midpoint, vers, tris, edges, EMAP, EF, EI, workQueue, EQ, C))   // collapse_edge()����2.2
            break; 
        FlagCollapsed = true; 
        num_collapsed++; 
      }

      if(FlagCollapsed)
      {
        viewer.data().clear(); 
        viewer.data().set_mesh(vers, tris); 
        viewer.data().set_face_based(true); 
      }
    }

    return false; 
  }; 

  // lambda���������¼���Ӧ��
  const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
  {
    switch(key)
    {
      case ' ':
        viewer.core().is_animating ^= 1; 
        break; 
      case 'R':
      case 'r':
        reset(); 
        break; 
      default:
        return false; 
    }
    return true; 
  }; 


  // 1. ��ʼ������λ����
  reset(); 

  // 2. visualize, working loop; 
  viewer.core().background_color.setConstant(1); 
  viewer.core().is_animating = true;                                 // ����
  viewer.callback_key_down = key_down; 
  viewer.callback_pre_draw = pre_draw; 
  return viewer.launch(); 
}
