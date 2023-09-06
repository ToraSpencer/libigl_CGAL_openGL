#ifndef IGL_ARAP_H
#define IGL_ARAP_H
#include "igl_inline.h"
#include "min_quad_with_fixed.h"
#include "ARAPEnergyType.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

// ARAP deformation(As Rigid As Possible Deformation)——尽可能刚性的变形；
namespace igl
{
  struct ARAPData
  {
      /*
         n  #vers
         G  #vers list of group indices (1 to k) for each vertex, such that vertex i is assigned to group G(i)
                    energy  type of energy to use  with_dynamics  whether using dynamics (need to call arap_precomputation
                    after changing)
         f_ext      #vers by dim list of external forces
         vel         #vers by dim list of velocities
         h            dynamics time step
         ym         ~Young's modulus smaller is softer, larger is more rigid/stiff
         max_iter   maximum inner iterations
         K          rhs pre-multiplier
         M          mass matrix
         solver_data    quadratic solver data
         b              list of boundary indices into vers
         dim        dimension being used for solving
      */
    int n;
    Eigen::VectorXi G;
    ARAPEnergyType energy;
    bool with_dynamics;
    Eigen::MatrixXd f_ext,vel;
    double h;
    double ym;
    int max_iter;
    Eigen::SparseMatrix<double> K,M;
    Eigen::SparseMatrix<double> CSM;
    min_quad_with_fixed_data<double> solver_data;
    Eigen::VectorXi b;
    int dim;
      ARAPData():
        n(0),
        G(),
        energy(ARAP_ENERGY_TYPE_DEFAULT),
        with_dynamics(false),
        f_ext(),
        h(1),
        ym(1),
        max_iter(10),
        K(),
        CSM(),
        solver_data(),
        b(),
        dim(-1) // force this to be set by _precomputation
    {
    };
  };
  

  // ARAP变形的预处理计算：
  /*
   Compute necessary information to start using an ARAP deformation
  
   Inputs:
     vers      #vers by dim list of mesh positions
     tris       #tris by simplex-size list of triangle|tet indices into vers
     dim      dimension being used at solve time. For deformation usually dim =
                            vers.cols(), for surface parameterization vers.cols() = 3 and dim = 2
     b          #b list of "boundary" fixed vertex indices into vers

   Outputs:
     data       struct containing necessary precomputation
  */
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedb>
  IGL_INLINE bool arap_precomputation(
    const Eigen::PlainObjectBase<DerivedV> & vers,
    const Eigen::PlainObjectBase<DerivedF> & tris,
    const int dim,
    const Eigen::PlainObjectBase<Derivedb> & b,
    ARAPData & data);


  // ARAP变形：
  /*
       Inputs:
         bc             输入的边界条件；#b by dim list of boundary conditions
         data          arap_precomputation()计算出的预处理数据；struct containing necessary precomputation and parameters
         versOut              输出点云；#vers by dim initial guess
  */
  template <
    typename Derivedbc,
    typename DerivedU>
  IGL_INLINE bool arap_solve(
    const Eigen::PlainObjectBase<Derivedbc> & bc,
    ARAPData& data,
    Eigen::PlainObjectBase<DerivedU>& versOut);
};

#ifndef IGL_STATIC_LIBRARY
#include "arap.cpp"
#endif

#endif
