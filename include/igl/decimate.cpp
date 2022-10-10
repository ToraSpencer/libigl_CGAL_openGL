#include "decimate.h"
#include "collapse_edge.h"
#include "edge_flaps.h"
#include "decimate_trivial_callbacks.h"
#include "is_edge_manifold.h"
#include "remove_unreferenced.h"
#include "slice_mask.h"
#include "slice.h"
#include "connect_boundary_to_infinity.h"
#include "parallel_for.h"
#include "max_faces_stopping_condition.h"
#include "shortest_edge_and_midpoint.h"


// ����1.1����ʹ��Ĭ�ϵĺ�����ִ�б��۵����񾫼�
IGL_INLINE bool igl::decimate(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const size_t max_m, 
  Eigen::MatrixXd & versOut, 
  Eigen::MatrixXi & trisOut, 
  Eigen::VectorXi & newOldTrisInfo, 
  Eigen::VectorXi & newOldVersInfo)
{
  const int trisCountOri = tris.rows();          // Original number of faces
  int trisCount = tris.rows();                          // Tracking number of faces

  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;

  igl::connect_boundary_to_infinity(vers, tris, VO, FO);
  Eigen::VectorXi edgeUeInfo;
  Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
  edge_flaps(FO, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);


  // note:
  /*
       decimate will not work correctly on non-edge-manifold meshes. 

       By extension this includes meshes with non-manifold vertices on the boundary since these
            will create a non-manifold edge when connected to infinity.
  */

  {
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
    Eigen::Array<bool, Eigen::Dynamic, 1> BE;
    if(!is_edge_manifold(FO, uEdges.rows(), edgeUeInfo, BF, BE))
        return false;
  }


  decimate_pre_collapse_callback always_try;
  decimate_post_collapse_callback never_care;
  decimate_trivial_callbacks(always_try, never_care);             // ����Ĭ�ϵ�pre_collapse��post_collapse�����ӣ�

  bool ret = decimate(VO, FO, 
        shortest_edge_and_midpoint,                                                                   // cost_and_placement()�����ӣ�                                                                 
        max_faces_stopping_condition(trisCount, trisCountOri, max_m),             // ����stopping_condition()�����ӣ�
        always_try, never_care,                                                             // pre_collapse��post_collapse������
        uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, 
        versOut, trisOut, 
        newOldTrisInfo, 
        newOldVersInfo);

  const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (newOldTrisInfo.array() < trisCountOri);
  igl::slice_mask(Eigen::MatrixXi(trisOut), keep, 1, trisOut);
  igl::slice_mask(Eigen::VectorXi(newOldTrisInfo), keep, 1, newOldTrisInfo);

  Eigen::VectorXi _1, I2;
  igl::remove_unreferenced(Eigen::MatrixXd(versOut), Eigen::MatrixXi(trisOut), versOut, trisOut, _1, I2);
  igl::slice(Eigen::VectorXi(newOldVersInfo), I2, 1, newOldVersInfo);

  return ret;
}


// ����1.1.1
IGL_INLINE bool igl::decimate(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const size_t max_m, 
  Eigen::MatrixXd & versOut, 
  Eigen::MatrixXi & trisOut, 
  Eigen::VectorXi & newOldTrisInfo)
{
  Eigen::VectorXi newOldVersInfo;
  return igl::decimate(vers, tris, max_m, versOut, trisOut, newOldTrisInfo, newOldVersInfo);
}


// ����1.2
IGL_INLINE bool igl::decimate(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_stopping_condition_callback & stopping_condition, 
  Eigen::MatrixXd & versOut, 
  Eigen::MatrixXi & trisOut, 
  Eigen::VectorXi & newOldTrisInfo, 
  Eigen::VectorXi & newOldVersInfo
  )
{
  decimate_pre_collapse_callback always_try;
  decimate_post_collapse_callback never_care;
  decimate_trivial_callbacks(always_try, never_care);
  return igl::decimate(
    vers, tris, cost_and_placement, stopping_condition, always_try, never_care, versOut, trisOut, newOldTrisInfo, newOldVersInfo);
}


// ����1.3
IGL_INLINE bool igl::decimate(
  const Eigen::MatrixXd & vers, 
  const Eigen::MatrixXi & tris, 
  const decimate_cost_and_placement_callback & cost_and_placement, 
  const decimate_stopping_condition_callback & stopping_condition, 
  const decimate_pre_collapse_callback       & pre_collapse, 
  const decimate_post_collapse_callback      & post_collapse, 
  Eigen::MatrixXd & versOut, 
  Eigen::MatrixXi & trisOut, 
  Eigen::VectorXi & newOldTrisInfo, 
  Eigen::VectorXi & newOldVersInfo
  )
{
  Eigen::VectorXi edgeUeInfo;
  Eigen::MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
  edge_flaps(tris, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);
  return igl::decimate(
    vers, tris, 
    cost_and_placement, stopping_condition, pre_collapse, post_collapse, 
    uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, 
    versOut, trisOut, newOldTrisInfo, newOldVersInfo);
}


// ����1.
IGL_INLINE bool igl::decimate(
    const Eigen::MatrixXd& vers,
    const Eigen::MatrixXi& tris,
    const decimate_cost_and_placement_callback& cost_and_placement,
    const decimate_stopping_condition_callback& stopping_condition,
    const decimate_pre_collapse_callback& pre_collapse,
    const decimate_post_collapse_callback& post_collapse,
    const Eigen::MatrixXi& OE,
    const Eigen::VectorXi& OEMAP,
    const Eigen::MatrixXi& OEF,
    const Eigen::MatrixXi& OEI,
    Eigen::MatrixXd& versOut,
    Eigen::MatrixXi& trisOut,
    Eigen::VectorXi& newOldTrisInfo,
    Eigen::VectorXi& newOldVersInfo)
{
    using namespace Eigen;
    using namespace std;

    // 1. ׼������
    Eigen::MatrixXd versCopy = vers;
    Eigen::MatrixXi trisCopy = tris;
    VectorXi edgeUeInfo;
    MatrixXi uEdges, UeTrisInfo, UeCornersInfo;
    edge_flaps(trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo);

    {
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
        Eigen::Array<bool, Eigen::Dynamic, 1> BE;
        if (!is_edge_manifold(trisCopy, uEdges.rows(), edgeUeInfo, BF, BE))
            return false;
    }

            // ���ȶ��У�
    igl::min_heap<std::tuple<double, int, int> > pQueue;

            // Could reserve with https://stackoverflow.com/a/29236236/148668
    Eigen::VectorXi EQ = Eigen::VectorXi::Zero(uEdges.rows());

            // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd collapsed(uEdges.rows(), versCopy.cols());

            // note
            /*
                 Pushing into a vector then using constructor was slower. 
                        Maybe using std::move + make_heap would squeeze out something?

                 Separating the cost/placement evaluation from the pQueue filling is a  performance hit for serial 
                        but faster if we can parallelize the cost/placement.
            */


    // 2. ����ÿ������ߵ�costֵ���Դ�Ϊ���ȼ��������ȶ���
    {
        Eigen::VectorXd costs(uEdges.rows());
        igl::parallel_for(uEdges.rows(), [&](const int e)
            {
                double cost = e;
                RowVectorXd p(1, 3);
                cost_and_placement(e, versCopy, trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, cost, p);
                collapsed.row(e) = p;
                costs(e) = cost;
            },
            10000);

        for (int i = 0; i < uEdges.rows(); i++)
            pQueue.emplace(costs(i), i, 0);
    }


    // 3. ���ȶ�����ִ�б��۵���
    int prev_e = -1;
    bool clean_finish = false;
    while (true)
    {
        int e, e1, e2, f1, f2;
        if (collapse_edge(cost_and_placement, pre_collapse, post_collapse,\
                versCopy, trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, \
                pQueue, EQ, collapsed, e, e1, e2, f1, f2))          // collapse_edge()����2��
        {
            if (stopping_condition(versCopy, trisCopy, uEdges, edgeUeInfo, UeTrisInfo, UeCornersInfo, pQueue, EQ, collapsed, e, e1, e2, f1, f2))
            {
                clean_finish = true;
                break;
            }
        }
        else
        {
            if (e == -1)
                break;                // a candidate edge was not even found in pQueue.
            if (prev_e == e)
            {
                assert(false && "Edge collapse no progress... bad stopping condition?");
                break;
            }

            // Edge was not collapsed... must have been invalid. collapse_edge should have updated its cost to inf... continue
        }
        prev_e = e;
    }


    // 4. ɾ�����к��б��ΪIGL_COLLAPSE_EDGE_NULL�ߵ�����Ƭ��
    MatrixXi tris0(trisCopy.rows(), 3);
    VectorXi _1;
    newOldTrisInfo.resize(trisCopy.rows());
    int m = 0;
    for (int i = 0; i < trisCopy.rows(); i++)
    {
        if (trisCopy(i, 0) != IGL_COLLAPSE_EDGE_NULL ||
            trisCopy(i, 1) != IGL_COLLAPSE_EDGE_NULL ||
            trisCopy(i, 2) != IGL_COLLAPSE_EDGE_NULL)
        {
            tris0.row(m) = trisCopy.row(i);
            newOldTrisInfo(m) = i;
            m++;
        }
    }
    tris0.conservativeResize(m, tris0.cols());              // �����൱��shrink_to_fit();
    newOldTrisInfo.conservativeResize(m);

    // 5. ɾ�������еĹ������㣺
    igl::remove_unreferenced(versCopy, tris0, versOut, trisOut, _1, newOldVersInfo);            // ����1.1

    return clean_finish;
}
