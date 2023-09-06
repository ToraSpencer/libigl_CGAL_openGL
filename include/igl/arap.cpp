// This file is part of libigl,  a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file,  You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "arap.h"
#include "colon.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "group_sum_matrix.h"
#include "covariance_scatter_matrix.h"
#include "speye.h"
#include "mode.h"
#include "project_isometrically_to_plane.h"
#include "slice.h"
#include "arap_rhs.h"
#include "repdiag.h"
#include "columnize.h"
#include "fit_rotations.h"
#include <cassert>
#include <iostream>

template <
  typename DerivedV, 
  typename DerivedF, 
  typename Derivedb>
IGL_INLINE bool igl::arap_precomputation(
  const Eigen::PlainObjectBase<DerivedV> & vers, 
  const Eigen::PlainObjectBase<DerivedF> & tris, 
  const int dim, 
  const Eigen::PlainObjectBase<Derivedb> & b, 
  ARAPData & data)
{
  using namespace std; 
  using namespace Eigen; 
  typedef typename DerivedV::Scalar Scalar;  

  const int n = vers.rows(); 
  data.n = n; 
  assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds"); 
  assert((b.size() == 0 || b.minCoeff() >=0) && "b out of bounds");  
  data.b = b;  
  assert((dim == 3 || dim ==2) && "dim should be 2 or 3"); 
  data.dim = dim;  
  // Defaults
  data.f_ext = MatrixXd::Zero(n, data.dim); 

  assert(data.dim <= vers.cols() && "solve dim should be <= embedding"); 
  bool flat = (vers.cols() - data.dim)==1; 

  DerivedV plane_V; 
  DerivedF plane_F; 
  typedef SparseMatrix<Scalar> SparseMatrixS; 
  SparseMatrixS ref_map, ref_map_dim; 
  if(flat)
  {
    project_isometrically_to_plane(vers, tris, plane_V, plane_F, ref_map); 
    repdiag(ref_map, dim, ref_map_dim); 
  }
  const PlainObjectBase<DerivedV>& ref_V = (flat?plane_V:vers); 
  const PlainObjectBase<DerivedF>& ref_F = (flat?plane_F:tris); 
  SparseMatrixS L; 
  cotmatrix(vers, tris, L); 

  ARAPEnergyType eff_energy = data.energy; 
  if(eff_energy == ARAP_ENERGY_TYPE_DEFAULT)
  {
    switch(tris.cols())
    {
      case 3:
        if(data.dim == 3) 
          eff_energy = ARAP_ENERGY_TYPE_SPOKES_AND_RIMS; 
        else 
          eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;  
        break; 
      case 4:
        eff_energy = ARAP_ENERGY_TYPE_ELEMENTS; 
        break; 
      default:
        assert(false); 
    }
  }

  // Get covariance scatter matrix,  when applied collects the covariance matrices used to fit rotations to during optimization
  covariance_scatter_matrix(ref_V, ref_F, eff_energy, data.CSM); 
  if(flat) 
    data.CSM = (data.CSM * ref_map_dim.transpose()).eval();  
  assert(data.CSM.cols() == vers.rows()*data.dim); 

  // Get group sum scatter matrix,  when applied sums all entries of the same group according to G
  SparseMatrix<double> G_sum; 
  if(data.G.size() == 0)
  {
    if(eff_energy == ARAP_ENERGY_TYPE_ELEMENTS) 
      speye(tris.rows(), G_sum); 
    else 
      speye(n, G_sum);  
  }
  else
  {
    // groups are defined per vertex,  convert to per face using mode
    if(eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
    {
      Eigen::Matrix<int, Eigen::Dynamic, 1> GG; 
      MatrixXi GF(tris.rows(), tris.cols()); 
      for(int j = 0; j<tris.cols(); j++)
      {
        Matrix<int, Eigen::Dynamic, 1> GFj; 
        slice(data.G, tris.col(j), GFj); 
        GF.col(j) = GFj; 
      }
      mode<int>(GF, 2, GG); 
      data.G=GG; 
    }
    //printf("group_sum_matrix()\n"); 
    group_sum_matrix(data.G, G_sum); 
  }
  SparseMatrix<double> G_sum_dim; 
  repdiag(G_sum, data.dim, G_sum_dim); 
  assert(G_sum_dim.cols() == data.CSM.rows()); 
  data.CSM = (G_sum_dim * data.CSM).eval(); 

  arap_rhs(ref_V, ref_F, data.dim, eff_energy, data.K); 
  if(flat) 
    data.K = (ref_map_dim * data.K).eval();  

  assert(data.K.rows() == data.n*data.dim); 

  SparseMatrix<double> Q = (-L).eval(); 

  if(data.with_dynamics)
  {
    const double h = data.h; 
    assert(h != 0); 
    SparseMatrix<double> M; 
    massmatrix(vers, tris, MASSMATRIX_TYPE_DEFAULT, data.M); 
    const double dw = (1./data.ym)*(h*h); 
    SparseMatrix<double> DQ = dw * 1./(h*h)*data.M; 
    Q += DQ; 
    // Dummy external forces
    data.f_ext = MatrixXd::Zero(n, data.dim); 
    data.vel = MatrixXd::Zero(n, data.dim); 
  }

  return min_quad_with_fixed_precompute(Q, b, SparseMatrix<double>(), true, data.solver_data); 
}


template <
  typename Derivedbc, 
  typename DerivedU>
IGL_INLINE bool igl::arap_solve(
  const Eigen::PlainObjectBase<Derivedbc> & bc, 
  ARAPData & data, 
  Eigen::PlainObjectBase<DerivedU> & versOut)
{
  using namespace Eigen; 
  using namespace std; 
  assert(data.b.size() == bc.rows()); 
  if(bc.size() > 0) 
    assert(bc.cols() == data.dim && "bc.cols() match data.dim");  

  const int n = data.n; 
  int iter = 0; 
  if(versOut.size() == 0)
  {
    // terrible initial guess.. should at least copy input mesh
#ifndef NDEBUG
    cerr<<"arap_solve: Using terrible initial guess for versOut. Try versOut = vers."<<endl; 
#endif
    versOut = MatrixXd::Zero(data.n, data.dim); 
  }
  else 
    assert(versOut.cols() == data.dim && "versOut.cols() match data.dim");  

  // changes each arap iteration
  MatrixXd versOut_prev = versOut; 

  // doesn't change for fixed with_dynamics timestep
  MatrixXd vers0; 
  if(data.with_dynamics) 
    vers0 = versOut_prev;  

  while(iter < data.max_iter)
  {
    versOut_prev = versOut; 
    // enforce boundary conditions exactly
    for(int bi = 0; bi<bc.rows(); bi++) 
      versOut.row(data.b(bi)) = bc.row(bi);  

    const auto & Udim = versOut.replicate(data.dim, 1); 
    assert(versOut.cols() == data.dim); 

    // As if versOut.col(2) was 0
    MatrixXd S = data.CSM * Udim; 
    // THIS NORMALIZATION IS IMPORTANT TO GET SINGLE PRECISION SVD CODE TO WORK CORRECTLY.
    S /= S.array().abs().maxCoeff(); 

    const int Rdim = data.dim; 
    MatrixXd R(Rdim, data.CSM.rows()); 
    if(R.rows() == 2) 
      fit_rotations_planar(S, R); 
    else
    {
      fit_rotations(S, true, R); 
//#ifdef __SSE__ // fit_rotations_SSE will convert to float if necessary
//      fit_rotations_SSE(S, R); 
//#else
//      fit_rotations(S, true, R); 
//#endif
    }  

    // Number of rotations: #vertices or #elements
    int num_rots = data.K.cols()/Rdim/Rdim; 

    // distribute group rotations to vertices in each group
    MatrixXd eff_R; 
    if(data.G.size() == 0) 
      eff_R = R; 
    else
    {
      eff_R.resize(Rdim, num_rots*Rdim); 
      for(int r = 0; r<num_rots; r++) 
        eff_R.block(0, Rdim*r, Rdim, Rdim) = R.block(0, Rdim*data.G(r), Rdim, Rdim);  
    }

    MatrixXd Dl; 
    if(data.with_dynamics)
    {
      assert(data.M.rows() == n && "No mass matrix. Call arap_precomputation if changing with_dynamics"); 
      const double h = data.h; 
      assert(h != 0); 
      //Dl = 1./(h*h*h)*M*(-2.*V0 + Vm1) - fext; 
      // data.vel = (V0-Vm1)/h
      // h*data.vel = (V0-Vm1)
      // -h*data.vel = -V0+Vm1)
      // -V0-h*data.vel = -2V0+Vm1
      const double dw = (1./data.ym)*(h*h); 
      Dl = dw * (1./(h*h)*data.M*(-vers0 - h*data.vel) - data.f_ext); 
    }

    VectorXd Rcol; 
    columnize(eff_R, num_rots, 2, Rcol); 
    VectorXd Bcol = -data.K * Rcol; 
    assert(Bcol.size() == data.n*data.dim); 
    for(int c = 0; c<data.dim; c++)
    {
      VectorXd Uc, Bc, bcc, Beq; 
      Bc = Bcol.block(c*n, 0, n, 1); 
      if(data.with_dynamics) 
        Bc += Dl.col(c);  
      if(bc.size()>0) 
        bcc = bc.col(c);  

      min_quad_with_fixed_solve(
        data.solver_data, 
        Bc, bcc, Beq, 
        Uc); 
      versOut.col(c) = Uc; 
    }

    iter++; 
  }

  if(data.with_dynamics) 
    data.vel = (versOut-vers0)/data.h;          // Keep track of velocity for next time 

  return true; 
}



#ifdef IGL_STATIC_LIBRARY
template bool igl::arap_solve<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  igl::ARAPData&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&); 
template bool igl::arap_precomputation<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  1,  0,  -1,  1> >(Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> > const&,  int,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  1,  0,  -1,  1> > const&,  igl::ARAPData&); 
#endif
