#include "simplify_polyhedron.h"
#include "decimate.h"
#include "circulation.h"
#include "per_face_normals.h"
#include "infinite_cost_stopping_condition.h"
#include <functional>

IGL_INLINE void igl::simplify_polyhedron(
  const Eigen::MatrixXd & OV,
  const Eigen::MatrixXi & OF,
  Eigen::MatrixXd & vers,
  Eigen::MatrixXi & tris,
  Eigen::VectorXi & J)
{
  // TODO: to generalize to open meshes, 0-cost should keep all incident
  // boundary edges on their original lines. (for non-manifold meshes,
  // igl::decimate needs to be generalized)

  Eigen::MatrixXd N;
  // Function for computing cost of collapsing edge (0 if at least one
  // direction doesn't change pointset, inf otherwise) and placement (in lowest
  // cost direction).
  const auto & perfect= [&N](
    const int e,
    const Eigen::MatrixXd & vers,
    const Eigen::MatrixXi & tris,
    const Eigen::MatrixXi & E,
    const Eigen::VectorXi & edgeUeInfo,
    const Eigen::MatrixXi & EF,
    const Eigen::MatrixXi & EI,
    double & cost,
    Eigen::RowVectorXd & p)
  {
    // Function for ocmputing cost (0 or inf) of collapsing edge by placing
    // vertex at `positive` end of edge.
    const auto & perfect_directed = [&N](
      const int e,
      const bool positive,
      const Eigen::MatrixXd & vers,
      const Eigen::MatrixXi & tris,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & edgeUeInfo,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      double & cost,
      Eigen::RowVectorXd & p)
    {
      const auto vi = E(e,positive);
      const auto vj = E(e,!positive);
      p = vers.row(vj);
      std::vector<int> faces = igl::circulation(e,positive,edgeUeInfo,EF,EI);
      cost = 0;
      for(auto f : faces)
      {
        // Skip the faces being collapsed
        if(f == EF(e,0) || f == EF(e,1))
            continue;
        const Eigen::RowVectorXd nbefore = N.row(f);
        // Face with vi replaced with vj
        const Eigen::RowVector3i fafter(
            tris(f,0) == vi ? vj : tris(f,0),
            tris(f,1) == vi ? vj : tris(f,1),
            tris(f,2) == vi ? vj : tris(f,2));
        Eigen::RowVectorXd nafter;
        igl::per_face_normals(vers,fafter,nafter);
        const double epsilon = 1e-10;
        // if normal changed then not feasible, break
        if((nbefore-nafter).norm() > epsilon)
        {
          cost = std::numeric_limits<double>::infinity();
          break;
        }
      }
    }; 
    p.resize(3);
    double cost0, cost1;
    Eigen::RowVectorXd p0, p1;
    perfect_directed(e,false,vers,tris,E,edgeUeInfo,EF,EI,cost0,p0);
    perfect_directed(e,true,vers,tris,E,edgeUeInfo,EF,EI,cost1,p1);
    if(cost0 < cost1)
    {
      cost = cost0;
      p = p0;
    }else
    {
      cost = cost1;
      p = p1;
    }
  };
  igl::per_face_normals(OV,OF,N);
  Eigen::VectorXi I;
  igl::decimate(
    OV,OF,
    perfect,
    igl::infinite_cost_stopping_condition(perfect),
    vers,tris,J,I);
}

