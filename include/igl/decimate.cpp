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


// 重载1.1
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
  int trisCount = tris.rows();                   // Tracking number of faces

  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;

  igl::connect_boundary_to_infinity(vers, tris, VO, FO);
  Eigen::VectorXi edgeUeInfo;
  Eigen::MatrixXi E, EF, EI;
  edge_flaps(FO, E, edgeUeInfo, EF, EI);

  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  {
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
    Eigen::Array<bool, Eigen::Dynamic, 1> BE;
    if(!is_edge_manifold(FO, E.rows(), edgeUeInfo, BF, BE))
        return false;
  }

  decimate_pre_collapse_callback always_try;
  decimate_post_collapse_callback never_care;
  decimate_trivial_callbacks(always_try, never_care);

  bool ret = decimate(VO, FO, 
        shortest_edge_and_midpoint, 
        max_faces_stopping_condition(trisCount, trisCountOri, max_m), 
        always_try, never_care, E, edgeUeInfo, EF, EI, 
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


// 重载1.1.1
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


// 重载1.2
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


// 重载1.3
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
  Eigen::MatrixXi E, EF, EI;
  edge_flaps(tris, E, edgeUeInfo, EF, EI);
  return igl::decimate(
    vers, tris, 
    cost_and_placement, stopping_condition, pre_collapse, post_collapse, 
    E, edgeUeInfo, EF, EI, 
    versOut, trisOut, newOldTrisInfo, newOldVersInfo);
}


// 重载1.
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

    // Working copies
    Eigen::MatrixXd versCopy = vers;
    Eigen::MatrixXi trisCopy = tris;
    VectorXi edgeUeInfo;
    MatrixXi E, EF, EI;
    edge_flaps(trisCopy, E, edgeUeInfo, EF, EI);

    {
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> BF;
        Eigen::Array<bool, Eigen::Dynamic, 1> BE;
        if (!is_edge_manifold(trisCopy, E.rows(), edgeUeInfo, BF, BE))
            return false;
    }

    igl::min_heap<std::tuple<double, int, int> > Q;

    // Could reserve with https://stackoverflow.com/a/29236236/148668
    Eigen::VectorXi EQ = Eigen::VectorXi::Zero(E.rows());

    // If an edge were collapsed, we'd collapse it to these points:
    MatrixXd C(E.rows(), versCopy.cols());

    // Pushing into a vector then using constructor was slower. Maybe using std::move + make_heap would squeeze out something?

    // Separating the cost/placement evaluation from the Q filling is a
    // performance hit for serial but faster if we can parallelize the cost/placement.
    {
        Eigen::VectorXd costs(E.rows());
        igl::parallel_for(E.rows(), [&](const int e)
            {
                double cost = e;
                RowVectorXd p(1, 3);
                cost_and_placement(e, versCopy, trisCopy, E, edgeUeInfo, EF, EI, cost, p);
                C.row(e) = p;
                costs(e) = cost;
            },
            10000);
        for (int e = 0; e < E.rows(); e++)
            Q.emplace(costs(e), e, 0);
    }

    int prev_e = -1;
    bool clean_finish = false;

    while (true)
    {
        int e, e1, e2, f1, f2;
        if (collapse_edge(cost_and_placement, pre_collapse, post_collapse,
            versCopy, trisCopy, E, edgeUeInfo, EF, EI, Q, EQ, C, e, e1, e2, f1, f2))
        {
            if (stopping_condition(versCopy, trisCopy, E, edgeUeInfo, EF, EI, Q, EQ, C, e, e1, e2, f1, f2))
            {
                clean_finish = true;
                break;
            }
        }
        else
        {
            if (e == -1)
                break;                // a candidate edge was not even found in Q.
            if (prev_e == e)
            {
                assert(false && "Edge collapse no progress... bad stopping condition?");
                break;
            }

            // Edge was not collapsed... must have been invalid. collapse_edge should have updated its cost to inf... continue
        }
        prev_e = e;
    }

    // remove all IGL_COLLAPSE_EDGE_NULL faces
    MatrixXi F2(trisCopy.rows(), 3);
    newOldTrisInfo.resize(trisCopy.rows());
    int m = 0;
    for (int f = 0; f < trisCopy.rows(); f++)
    {
        if (
            trisCopy(f, 0) != IGL_COLLAPSE_EDGE_NULL ||
            trisCopy(f, 1) != IGL_COLLAPSE_EDGE_NULL ||
            trisCopy(f, 2) != IGL_COLLAPSE_EDGE_NULL)
        {
            F2.row(m) = trisCopy.row(f);
            newOldTrisInfo(m) = f;
            m++;
        }
    }
    F2.conservativeResize(m, F2.cols());
    newOldTrisInfo.conservativeResize(m);
    VectorXi _1;
    igl::remove_unreferenced(versCopy, F2, versOut, trisOut, _1, newOldVersInfo);

    return clean_finish;
}
