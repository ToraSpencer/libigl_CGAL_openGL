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
  VectorXi edgeUeInfo;
  MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
  Eigen::VectorXi timeStamps;
  MatrixXd collapsedVers;                       // ������֮�����ɵĶ����λ�ã�
  int num_collapsed;                                // ������ѭ���ļ�����
  igl::min_heap< std::tuple<double, int, int>> pQueue;
  igl::opengl::glfw::Viewer viewer; 
  read_triangle_mesh(filename, versOld, trisOld); 


  // lambda�����������ݸ�λ
  const auto & reset = [&]()
  {
    tris = trisOld; 
    vers = versOld; 
    edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo); 
    collapsedVers.resize(uEdges.rows(), vers.cols()); 
    VectorXd costs(uEdges.rows()); 

    pQueue = {};            // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
    timeStamps = Eigen::VectorXi::Zero(uEdges.rows()); 

    // r1. ����ÿ���ߵ��۵�costֵ���Լ��۵�֮��Ķ������꣬�������ȶ���pQueue��
    {
      Eigen::VectorXd costs(uEdges.rows()); 
      igl::parallel_for(uEdges.rows(), \
          [&](const int i)
          {
            double cost = i; 
            RowVectorXd edgeCenter(1, 3);           // ȡ�ߵ��е���Ϊ���۵�֮��Ķ��㣻
            shortest_edge_and_midpoint(i, vers, tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, edgeCenter); 
            collapsedVers.row(i) = edgeCenter; 
            costs(i) = cost; 
          },\
          10000); 

      for(int i = 0; i < uEdges.rows(); i++)
        pQueue.emplace(costs(i), i, 0); 
    }

    num_collapsed = 0; 
    viewer.data().clear(); 
    viewer.data().set_mesh(vers, tris); 
    viewer.data().set_face_based(true); 
  }; 


  // lambda����ִ�����񾫼�׼����Ⱦ�����ݣ�
  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer)->bool
  {
    // p1. ÿһ�ζ���ѭ���У�����10%�ıߣ�
    if(viewer.core().is_animating && !pQueue.empty())
    {
      bool FlagCollapsed = false;           // ����ѭ���б������Ƿ�ִ�гɹ���

      // p1.1 ִ�б���������collapse edge
      const int max_iter = std::ceil(0.01*pQueue.size()); 
      for(int j = 0; j < max_iter; j++)
      {
        if(!collapse_edge(shortest_edge_and_midpoint, vers, tris, uEdges, \
                    edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, timeStamps, collapsedVers))              // collapse_edge()����2.2
            break; 
        FlagCollapsed = true; 
        num_collapsed++; 
      }

      // p1.2 
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
