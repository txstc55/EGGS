// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "slim.h"

#include "boundary_loop.h"
#include "cotmatrix.h"
#include "edge_lengths.h"
#include "grad.h"
#include "local_basis.h"
#include "repdiag.h"
#include "vector_area_matrix.h"
#include "arap.h"
#include "cat.h"
#include "doublearea.h"
#include "grad.h"
#include "local_basis.h"
#include "per_face_normals.h"
#include "slice_into.h"
#include "volume.h"
#include "polar_svd.h"
#include "flip_avoiding_line_search.h"
#include "mapping_energy_with_jacobians.h"

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#include "Timer.h"
#include "sparse_cached.h"
#include "AtA_cached.h"

#ifdef CHOLMOD
#include <Eigen/CholmodSupport>
#endif

// #undef SLIM_CACHED

#define USE_MKL_SOLVER

#include "inline_expansion/utils.hpp"

#include "ie_helper.hpp"

#include "test_record_util.h"
namespace igl
{
namespace slim
{
// Definitions of internal functions
IGL_INLINE void buildRhs(igl::SLIMData &s, const Eigen::SparseMatrix<double> &A);
IGL_INLINE void add_soft_constraints(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void add_soft_constraints(igl::SLIMData &s, Eigen::SparseMatrix<double, Eigen::RowMajor> &L);
IGL_INLINE void add_soft_constraints_numeric_together(igl::SLIMData &s);
IGL_INLINE void add_soft_constraints_numeric_seperate_pre(igl::SLIMData &s, Eigen::SparseMatrix<ie::NumericType> &L);
IGL_INLINE void add_soft_constraints_numeric_seperate_final(igl::SLIMData &s);
IGL_INLINE double compute_energy(igl::SLIMData &s, Eigen::MatrixXd &V_new);
IGL_INLINE double compute_soft_const_energy(igl::SLIMData &s,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F,
                                            Eigen::MatrixXd &V_o);

IGL_INLINE void solve_weighted_arap(igl::SLIMData &s,
                                    const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    Eigen::MatrixXd &uv,
                                    Eigen::VectorXi &soft_b_p,
                                    Eigen::MatrixXd &soft_bc_p);
IGL_INLINE void update_weights_and_closest_rotations(igl::SLIMData &s,
                                                     Eigen::MatrixXd &uv);
IGL_INLINE void compute_jacobians(igl::SLIMData &s, const Eigen::MatrixXd &uv);
IGL_INLINE void build_linear_system(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void build_linear_system_eigen(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void build_linear_system_cached(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void build_linear_system_mkl(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void build_linear_system_numeric_together(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void build_linear_system_numeric_seperate(igl::SLIMData &s, Eigen::SparseMatrix<double> &L);
IGL_INLINE void pre_calc(igl::SLIMData &s);

std::ofstream result_file;

// Implementation

IGL_INLINE void compute_jacobians(igl::SLIMData &s, const Eigen::MatrixXd &uv)
{
  if (s.F.cols() == 3)
  {
    // Ji=[D1*u,D2*u,D1*v,D2*v];
    s.Ji.col(0) = s.Dx * uv.col(0);
    s.Ji.col(1) = s.Dy * uv.col(0);
    s.Ji.col(2) = s.Dx * uv.col(1);
    s.Ji.col(3) = s.Dy * uv.col(1);
  }
  else /*tet mesh*/
  {
    // Ji=[D1*u,D2*u,D3*u, D1*v,D2*v, D3*v, D1*w,D2*w,D3*w];
    s.Ji.col(0) = s.Dx * uv.col(0);
    s.Ji.col(1) = s.Dy * uv.col(0);
    s.Ji.col(2) = s.Dz * uv.col(0);
    s.Ji.col(3) = s.Dx * uv.col(1);
    s.Ji.col(4) = s.Dy * uv.col(1);
    s.Ji.col(5) = s.Dz * uv.col(1);
    s.Ji.col(6) = s.Dx * uv.col(2);
    s.Ji.col(7) = s.Dy * uv.col(2);
    s.Ji.col(8) = s.Dz * uv.col(2);
  }
}

IGL_INLINE void update_weights_and_closest_rotations(igl::SLIMData &s, Eigen::MatrixXd &uv)
{
  compute_jacobians(s, uv);
  slim_update_weights_and_closest_rotations_with_jacobians(s.Ji, s.slim_energy, s.exp_factor, s.W, s.Ri);
}

IGL_INLINE void solve_weighted_arap(igl::SLIMData &s,
                                    const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    Eigen::MatrixXd &uv,
                                    Eigen::VectorXi &soft_b_p,
                                    Eigen::MatrixXd &soft_bc_p)
{
  if (s.first_called)
  {
    result_file.open("result.txt");
  }
  using namespace Eigen;

  Eigen::SparseMatrix<double> L;
  switch (s.method_type)
  {
  case 0:
    build_linear_system_eigen(s, L);
    break;
  case 1:
    build_linear_system_cached(s, L);
    break;
  case 2:
    build_linear_system_mkl(s, L);
    break;
  case 3:
    build_linear_system_numeric_together(s, L);
    break;
  case 4:
    build_linear_system_numeric_seperate(s, L);
    break;
  default:
    build_linear_system_eigen(s, L);
    break;
  }

  //t.start();
  // solve
  Eigen::VectorXd Uc;
  std::vector<double> x;
  if (s.result_vector.size() == 0)
  {
    x.resize(L.rows());
  }
  else
  {
    x.resize(s.n_rows);
  }
#ifndef CHOLMOD
  if (s.dim == 2)
  {
    igl::Timer t;
#ifndef USE_MKL_SOLVER
    SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    t.start();
    Uc = solver.compute(L).solve(s.rhs);
    t.stop();
    std::cout << "\n"
              << "Eigen solver took " << t.getElapsedTimeInMicroSec();
#else
    if (s.first_called)
    {
      t.start();
      pardiso_init(s.pardiso_data);
      if (s.result_vector.size() == 0)
      {
        pardiso_support_matrix(s.pardiso_data, L);
      }
      else
      {
        pardiso_support_matrix(s.pardiso_data, s.L_outer.data(), s.L_inner.data(), s.result_vector.data(), s.n_rows);
      }
      t.stop();
      std::cout << "\n"
                << "Init everything for the mkl solver took " << t.getElapsedTimeInMicroSec();
      t.start();
      pardiso_symbolic_factor(s.pardiso_data);
      t.stop();
      std::cout << "\n"
                << "Symbolic factorization took " << t.getElapsedTimeInMicroSec();
    }
    else
    {
      if (s.result_vector.size() == 0)
      {
        pardiso_support_value(s.pardiso_data, L.valuePtr());
      }
      else
      {
        pardiso_support_value(s.pardiso_data, s.result_vector.data());
      }
    }
    t.start();
    pardiso_numeric_factor(s.pardiso_data);
    t.stop();
    write_to_file(result_file, "Numeric factorization", t.getElapsedTimeInMicroSec(), s.first_called);
    std::vector<double> b(s.rhs.data(), s.rhs.data() + s.rhs.rows() * s.rhs.cols());
    t.start();
    pardiso_solve(s.pardiso_data, x.data(), b.data());
    t.stop();
    write_to_file(result_file, "MKL Pardiso solver", t.getElapsedTimeInMicroSec(), s.first_called);
    s.first_called = false;
#endif
  }
  else
  {
    igl::Timer t;
#ifndef USE_MKL_SOLVER
    // seems like CG performs much worse for 2D and way better for 3D
    Eigen::VectorXd guess(uv.rows() * s.dim);
    for (int i = 0; i < s.v_num; i++)
      for (int j = 0; j < s.dim; j++)
        guess(uv.rows() * j + i) = uv(i, j); // flatten vector

    ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper> cg;
    t.start();
    cg.setTolerance(1e-8);
    cg.compute(L);
    Uc = cg.solveWithGuess(s.rhs, guess);
    t.stop();
    std::cout << "\n"
              << "Eigen ConjugateGradient solver took " << t.getElapsedTimeInMicroSec();
#else
    if (s.first_called)
    {
      t.start();
      pardiso_init(s.pardiso_data);
      if (s.result_vector.size() == 0)
      {
        pardiso_support_matrix(s.pardiso_data, L);
      }
      else
      {
        pardiso_support_matrix(s.pardiso_data, s.L_outer.data(), s.L_inner.data(), s.result_vector.data(), s.n_rows);
      }
      t.stop();
      std::cout << "\nInit everything for the mkl solver took " << t.getElapsedTimeInMicroSec();
      t.start();
      pardiso_symbolic_factor(s.pardiso_data);
      t.stop();
      std::cout << "\nSymbolic factorization took " << t.getElapsedTimeInMicroSec();
    }
    else
    {
      if (s.result_vector.size() == 0)
      {
        pardiso_support_value(s.pardiso_data, L.valuePtr());
      }
      else
      {
        pardiso_support_value(s.pardiso_data, s.result_vector.data());
      }
    }
    t.start();
    pardiso_numeric_factor(s.pardiso_data);
    t.stop();
    write_to_file(result_file, "Numeric factorization", t.getElapsedTimeInMicroSec(), s.first_called);
    std::vector<double> b(s.rhs.data(), s.rhs.data() + s.rhs.rows() * s.rhs.cols());
    t.start();
    pardiso_solve(s.pardiso_data, x.data(), b.data());
    t.stop();
    write_to_file(result_file, "MKL Pardiso solver", t.getElapsedTimeInMicroSec(), s.first_called);
    s.first_called = false;
#endif
  }
#else
  CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  Uc = solver.compute(L).solve(s.rhs);
#endif
  for (int i = 0; i < s.dim; i++)
#ifndef USE_MKL_SOLVER
    uv.col(i) = Uc.block(i * s.v_n, 0, s.v_n, 1);
#else
    uv.col(i) = Eigen::Map<Eigen::VectorXd>(x.data() + i * s.v_n, s.v_n);
#endif

  // t.stop();
  // std::cerr << "solve: " << t.getElapsedTime() << std::endl;
} // namespace slim

IGL_INLINE void pre_calc(igl::SLIMData &s)
{
  if (!s.has_pre_calc)
  {
    s.v_n = s.v_num;
    s.f_n = s.f_num;

    if (s.F.cols() == 3)
    {
      s.dim = 2;
      Eigen::MatrixXd F1, F2, F3;
      igl::local_basis(s.V, s.F, F1, F2, F3);
      Eigen::SparseMatrix<double> G;
      igl::grad(s.V, s.F, G);
      Eigen::SparseMatrix<double> Face_Proj;

      auto face_proj = [](Eigen::MatrixXd &F) {
        std::vector<Eigen::Triplet<double>> IJV;
        int f_num = F.rows();
        for (int i = 0; i < F.rows(); i++)
        {
          IJV.push_back(Eigen::Triplet<double>(i, i, F(i, 0)));
          IJV.push_back(Eigen::Triplet<double>(i, i + f_num, F(i, 1)));
          IJV.push_back(Eigen::Triplet<double>(i, i + 2 * f_num, F(i, 2)));
        }
        Eigen::SparseMatrix<double> P(f_num, 3 * f_num);
        P.setFromTriplets(IJV.begin(), IJV.end());
        return P;
      };

      s.Dx = face_proj(F1) * G;
      s.Dy = face_proj(F2) * G;
    }
    else
    {
      s.dim = 3;
      Eigen::SparseMatrix<double> G;
      igl::grad(s.V, s.F, G,
                s.mesh_improvement_3d /*use normal gradient, or one from a "regular" tet*/);
      s.Dx = G.block(0, 0, s.F.rows(), s.V.rows());
      s.Dy = G.block(s.F.rows(), 0, s.F.rows(), s.V.rows());
      s.Dz = G.block(2 * s.F.rows(), 0, s.F.rows(), s.V.rows());
    }

    s.W.resize(s.f_n, s.dim * s.dim);
    s.Dx.makeCompressed();
    s.Dy.makeCompressed();
    s.Dz.makeCompressed();
    s.Ri.resize(s.f_n, s.dim * s.dim);
    s.Ji.resize(s.f_n, s.dim * s.dim);
    s.rhs.resize(s.dim * s.v_num);

    // flattened weight matrix
    s.WGL_M.resize(s.dim * s.dim * s.f_n);
    for (int i = 0; i < s.dim * s.dim; i++)
      for (int j = 0; j < s.f_n; j++)
        s.WGL_M(i * s.f_n + j) = s.M(j);

    s.first_solve = true;
    s.has_pre_calc = true;
  }
}

IGL_INLINE void build_linear_system_eigen(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  result_file << "START BUILDING LINEAR SYSTEM USING EIGEN\n";
  std::vector<Eigen::Triplet<double>> IJV;
  igl::Timer t;
  t.start();
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  t.stop();
  write_to_file(result_file, "SLIM_BUILDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> A(s.dim * s.dim * s.f_n, s.dim * s.v_n);
  A.setFromTriplets(IJV.begin(), IJV.end());
  A.makeCompressed();
  t.stop();
  write_to_file(result_file, "CONSTRUCTING MATRIX A", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> At = A.transpose();
  At.makeCompressed();
  t.stop();
  write_to_file(result_file, "TAKING A TRANSPOSE", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> id_m(A.cols(), A.cols());
  id_m.setIdentity();
  t.stop();
  write_to_file(result_file, "SETTING ID_M", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  L = At * s.WGL_M.asDiagonal() * A + s.proximal_p * id_m; //add also a proximal term
  L.makeCompressed();
  t.stop();
  write_to_file(result_file, "COMPUTE L WITHOUT SOFT CONSTRAINTS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  buildRhs(s, A);
  t.stop();
  write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  add_soft_constraints(s, L);
  t.stop();
  t.start();
  L = L.triangularView<Eigen::Lower>();
  t.stop();
  write_to_file(result_file, "GETTING TRIANGULAR VIEW TOOK ", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void build_linear_system_cached(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  result_file << "START BUILDING LINEAR SYSTEM USING DANIELE'S METHOD\n";
  std::vector<Eigen::Triplet<double>> IJV;
  igl::Timer t;
  t.start();
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  t.stop();
  write_to_file(result_file, "SLIM_BUILDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  if (s.A.rows() == 0)
  {
    s.A = Eigen::SparseMatrix<double>(s.dim * s.dim * s.f_n, s.dim * s.v_n);
    igl::sparse_cached_precompute(IJV, s.A_data, s.A);
    // export_mtx(s.A, "small_face.mtx");
  }
  else
    igl::sparse_cached(IJV, s.A_data, s.A);
  t.stop();
  write_to_file(result_file, "CONSTRUCTING MATRIX A", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> id_m(s.A.cols(), s.A.cols());
  id_m.setIdentity();
  t.stop();
  write_to_file(result_file, "SETTING ID_M", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  s.AtA_data.W = s.WGL_M;
  if (s.AtA.rows() == 0)
  {
    igl::AtA_cached_precompute(s.A, s.AtA_data, s.AtA);
  }
  else
  {
    igl::AtA_cached(s.A, s.AtA_data, s.AtA);
  }
  t.stop();
  write_to_file(result_file, "COMPUTE ATDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  L = s.AtA + s.proximal_p * id_m; //add also a proximal
  L.makeCompressed();
  t.stop();
  write_to_file(result_file, "ADDING PROXIMAL TOOK ", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  buildRhs(s, s.A);
  t.stop();
  write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  add_soft_constraints(s, L);
  t.stop();
  t.start();
  L = L.triangularView<Eigen::Lower>();
  t.stop();
  write_to_file(result_file, "GETTING TRIANGULAR VIEW TOOK ", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void build_linear_system_mkl(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  result_file << "START BUILDING LINEAR SYSTEM USING MKL\n";
  std::vector<Eigen::Triplet<double>> IJV;
  igl::Timer t;
  t.start();
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  t.stop();
  write_to_file(result_file, "SLIM_BUILDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  if (s.A.rows() == 0)
  {
    s.A = Eigen::SparseMatrix<double>(s.dim * s.dim * s.f_n, s.dim * s.v_n);
    igl::sparse_cached_precompute(IJV, s.A_data, s.A);
    // export_mtx(s.A, "small_face.mtx");
  }
  else
    igl::sparse_cached(IJV, s.A_data, s.A);
  t.stop();
  write_to_file(result_file, "CONSTRUCTING MATRIX A", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  buildRhs(s, s.A);
  t.stop();
  write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  sparse_matrix_t A_mkl;
  create_mkl_csr_matrix(s.A, &A_mkl);
  t.stop();
  write_to_file(result_file, "CONVERTING A TO MKL FORMAT", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> wgl_m(s.A.rows(), s.A.rows());
  std::vector<Eigen::Triplet<double>> wgl_m_trip(s.A.rows());
  for (int i = 0; i < s.A.rows(); i++)
  {
    wgl_m_trip[i] = Eigen::Triplet<double>(i, i, s.WGL_M(i));
  }
  wgl_m.setFromTriplets(wgl_m_trip.begin(), wgl_m_trip.end());
  wgl_m.makeCompressed();
  sparse_matrix_t wgl_m_mkl;
  create_mkl_csr_matrix(wgl_m, &wgl_m_mkl);
  t.stop();
  write_to_file(result_file, "CREATING WGL_M_MKL TOOK ", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  Eigen::SparseMatrix<double> id_m(s.A.cols(), s.A.cols());
  id_m.setIdentity();
  id_m *= s.proximal_p;
  t.stop();
  write_to_file(result_file, "SETTING ID_M", t.getElapsedTimeInMicroSec(), s.first_called);
  add_soft_constraints(s, id_m);
  id_m.makeCompressed();
  t.start();
  sparse_matrix_t id_m_mkl;
  create_mkl_csr_matrix(id_m, &id_m_mkl);
  write_to_file(result_file, "CREATING PROXIMAL+SOFT_CONSTRAINTS MKL", t.getElapsedTimeInMicroSec(), s.first_called);
  t.stop();
  if (s.first_called)
  {
    mkl_sypr_pre(A_mkl, wgl_m_mkl, s.symmetric_type, &s.ata_mkl);
  }
  t.start();
  mkl_sypr_final(A_mkl, wgl_m_mkl, s.symmetric_type, &s.ata_mkl);
  t.stop();
  write_to_file(result_file, "COMPUTE ATDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  mkl_sparse_add(s.ata_mkl, id_m_mkl, 1, &s.final_result);
  t.stop();
  write_to_file(result_file, "ADDING PROXIMAL AND SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  if (s.first_called)
  {
    s.n_rows = s.A.cols();
    export_csr_from_mkl(s.final_result, s.L_outer, s.L_inner, s.result_vector);
  }
  else
  {
    export_val_from_mkl(s.final_result, s.result_vector);
  }
  t.stop();
  write_to_file(result_file, "EXPORTING CSR TOOK ", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void build_linear_system_numeric_together(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  result_file << "START BUILDING LINEAR SYSTEM USING NUMERIC DOING ALL COMPUTATIONS TOGETHER\n";
  std::vector<Eigen::Triplet<double>> IJV;
  igl::Timer t;
  t.start();
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  t.stop();
  write_to_file(result_file, "SLIM_BUILDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  if (s.A.rows() == 0)
  {
    s.A = Eigen::SparseMatrix<double>(s.dim * s.dim * s.f_n, s.dim * s.v_n);
    igl::sparse_cached_precompute(IJV, s.A_data, s.A);
    // export_mtx(s.A, "small_face.mtx");
  }
  else
    igl::sparse_cached(IJV, s.A_data, s.A);
  t.stop();
  write_to_file(result_file, "CONSTRUCTING MATRIX A", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  buildRhs(s, s.A);
  t.stop();
  write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  add_soft_constraints_numeric_together(s);
  t.stop();
  if (s.first_called)
  {
    s.n_rows = s.A.cols();
    // build A in numeric type
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> A_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(s.A, 0);
    s.datas[0].resize(s.A.nonZeros());
    s.datas[2].resize(s.WGL_M.rows() * s.WGL_M.cols());

    // build the soft constraints
    Eigen::SparseMatrix<double, Eigen::ColMajor> soft_constraints(s.A.cols(), s.A.cols());
    soft_constraints.setFromTriplets(s.soft_constraints_triplet.begin(), s.soft_constraints_triplet.end());
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> soft_constraints_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(soft_constraints, 1);

    // build the two diagonal matrices
    std::vector<Eigen::Triplet<ie::NumericType>> diagonals;
    diagonals.reserve(s.A.cols());
    ie::NumericType proximal = ie::NumericType(s.proximal_p);
    Eigen::SparseMatrix<double, Eigen::ColMajor> wgl_m(s.A.rows(), s.A.rows());
    std::vector<Eigen::Triplet<double>> wgl_m_trip;
    wgl_m_trip.reserve(s.A.rows());
    for (int i = 0; i < s.A.cols(); i++)
    {
      diagonals.push_back(Eigen::Triplet<ie::NumericType>(i, i, proximal));
    }
    for (int i = 0; i < s.A.rows(); i++)
    {
      wgl_m_trip.push_back(Eigen::Triplet<double>(i, i, 2.0));
    }
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> diagonal_numeric(s.A.cols(), s.A.cols());
    diagonal_numeric.setFromTriplets(diagonals.begin(), diagonals.end());
    wgl_m.setFromTriplets(wgl_m_trip.begin(), wgl_m_trip.end());
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> wgl_m_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(wgl_m, 2);

    // the actual structure initialization
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> result_numeric = Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor>(A_numeric.transpose()) * wgl_m_numeric * A_numeric + diagonal_numeric + soft_constraints_numeric;

    result_numeric = result_numeric.triangularView<Eigen::Lower>();
    s.ex = ie::NumericExecutor(result_numeric, 0);
    s.result_vector.resize(result_numeric.nonZeros(), 0);

    // copy the outer index pointer and inner index pointer to s
    s.L_outer.resize(result_numeric.rows() + 1);
    s.L_inner.resize(result_numeric.nonZeros());
    s.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
    s.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
    s.soft_constraints_triplet.resize(0);
    result_numeric.resize(0, 0);
  }
  t.start();
  s.datas[0].assign(s.A.valuePtr(), s.A.valuePtr() + s.A.nonZeros());
  s.datas[2].assign(s.WGL_M.data(), s.WGL_M.data() + s.WGL_M.rows() * s.WGL_M.cols());
  t.stop();
  write_to_file(result_file, "ASSIGNING DATAS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  s.ex.ExecuteMulti(s.datas, s.result_vector);
  t.stop();
  write_to_file(result_file, "COMPUTE EVERYTHING", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void build_linear_system_numeric_seperate(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  result_file << "START BUILDING LINEAR SYSTEM USING NUMERIC DOING ALL COMPUTATIONS SEPERATE\n";
  std::vector<Eigen::Triplet<double>> IJV;
  igl::Timer t;
  t.start();
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  t.stop();
  write_to_file(result_file, "SLIM_BUILDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  if (s.A.rows() == 0)
  {
    s.A = Eigen::SparseMatrix<double>(s.dim * s.dim * s.f_n, s.dim * s.v_n);
    igl::sparse_cached_precompute(IJV, s.A_data, s.A);
  }
  else
    igl::sparse_cached(IJV, s.A_data, s.A);
  t.stop();
  write_to_file(result_file, "CONSTRUCTING MATRIX A", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  buildRhs(s, s.A);
  t.stop();
  write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), s.first_called);
  if (s.first_called)
  {
    s.n_rows = s.A.cols();
    s.datas.resize(2);
    // build A in numeric type
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> A_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(s.A, 0);
    s.datas[0].resize(s.A.nonZeros());
    s.datas[1].resize(s.WGL_M.rows() * s.WGL_M.cols());

    // build wgl_m
    Eigen::SparseMatrix<double, Eigen::ColMajor> wgl_m(s.A.rows(), s.A.rows());
    std::vector<Eigen::Triplet<double>> wgl_m_trip;
    wgl_m_trip.reserve(s.A.rows());
    for (int i = 0; i < s.A.rows(); i++)
    {
      wgl_m_trip.push_back(Eigen::Triplet<double>(i, i, 2.0));
    }
    wgl_m.setFromTriplets(wgl_m_trip.begin(), wgl_m_trip.end());
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> wgl_m_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(wgl_m, 1);

    // the actual structure initialization
    Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> result_numeric = Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor>(A_numeric.transpose()) * wgl_m_numeric * A_numeric;
    result_numeric = result_numeric.triangularView<Eigen::Lower>();
    std::vector<int> zero_diagonal_row;
    // get the 0 diagonals
    for (int i = 0; i < s.A.cols(); i++)
    {
      if (s.A.outerIndexPtr()[i + 1] - s.A.outerIndexPtr()[i] == 0)
      {
        zero_diagonal_row.push_back(i);
      }
    }
    // adding proximal terms to the diagonals
    for (int i = 0; i < zero_diagonal_row.size(); i++)
    {
      result_numeric.insert(zero_diagonal_row[i], zero_diagonal_row[i]) = ie::NumericType(s.proximal_p);
    }

    result_numeric.makeCompressed();
    add_soft_constraints_numeric_seperate_pre(s, result_numeric);
    s.ex = ie::NumericExecutor(result_numeric, 0);
    s.result_vector.resize(result_numeric.nonZeros(), 0);

    // copy the outer index pointer and inner index pointer to s
    s.L_outer.resize(result_numeric.rows() + 1);
    s.L_inner.resize(result_numeric.nonZeros());
    s.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
    s.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
    result_numeric.resize(0, 0);
  }
  t.start();
  s.datas[0].assign(s.A.valuePtr(), s.A.valuePtr() + s.A.nonZeros());
  s.datas[1].assign(s.WGL_M.data(), s.WGL_M.data() + s.WGL_M.rows() * s.WGL_M.cols());
  t.stop();
  write_to_file(result_file, "ASSIGNING DATAS", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  s.ex.ExecuteMulti(s.datas, s.result_vector);
  t.stop();
  write_to_file(result_file, "COMPUTE ATDA", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  add_soft_constraints_numeric_seperate_final(s);
  t.stop();
}

IGL_INLINE void add_soft_constraints(igl::SLIMData &s, Eigen::SparseMatrix<double> &L)
{
  igl::Timer t;
  t.start();
  int v_n = s.v_num;
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
    }
  }
  t.stop();
  write_to_file(result_file, "RHS ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      L.coeffRef(d * v_n + v_idx, d * v_n + v_idx) += s.soft_const_p; // diagonal of matrix
    }
  }
  write_to_file(result_file, "L ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void add_soft_constraints(igl::SLIMData &s, Eigen::SparseMatrix<double, Eigen::RowMajor> &L)
{
  igl::Timer t;
  t.start();
  int v_n = s.v_num;
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
    }
  }
  t.stop();
  write_to_file(result_file, "RHS ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      L.coeffRef(d * v_n + v_idx, d * v_n + v_idx) += s.soft_const_p; // diagonal of matrix
    }
  }
  t.stop();
  write_to_file(result_file, "L ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void add_soft_constraints_numeric_together(igl::SLIMData &s)
{
  int v_n = s.v_num;
  int pos = 0;
  if (s.first_called)
  {
    s.datas.resize(3);
    std::vector<std::pair<int, int>> positions;
    positions.reserve(s.dim * s.b.rows());
    s.soft_constraints_triplet.reserve(s.dim * s.b.rows());
    s.datas[1].resize(s.dim * s.b.rows());

    for (int d = 0; d < s.dim; d++)
    {
      for (int i = 0; i < s.b.rows(); i++)
      {
        int v_idx = s.b(i);
        positions.push_back({d * v_n + v_idx, pos});
        s.soft_constraints_triplet.push_back(Eigen::Triplet<double>(d * v_n + v_idx, d * v_n + v_idx, 2.0));
        pos++;
      }
    }
    std::sort(positions.begin(), positions.end()); // sort the positions, because we want it to be the value array in csr format
    s.soft_constraints_pos.resize(positions.size());
    for (int i = 0; i < positions.size(); i++)
    {
      s.soft_constraints_pos[positions[i].second] = i;
    }
  }
  // now rearrange the values
  igl::Timer t;
  t.start();
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
    }
  }
  t.stop();
  write_to_file(result_file, "RHS ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  pos = 0;
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      s.datas[1][s.soft_constraints_pos[pos]] = s.soft_const_p;
      pos++;
    }
  }
  t.stop();
  write_to_file(result_file, "NUMERIC SET DATA 2", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE void add_soft_constraints_numeric_seperate_pre(igl::SLIMData &s, Eigen::SparseMatrix<ie::NumericType> &L)
{
  int v_n = s.v_num;
  std::vector<int> rows;
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      rows.push_back(d * v_n + v_idx);
    }
  }
  s.soft_constraints_index.resize(rows.size());
  int index = 0;
  for (int i = 0; i < L.outerSize(); ++i)
  {
    for (Eigen::SparseMatrix<ie::NumericType>::InnerIterator it(L, i); it; ++it)
    {
      if (it.row() == it.col())
      {
        if (it.value().operation != 1)
        {
          s.diagonal_index.push_back(index);
        }
        for (int j = 0; j < rows.size(); j++)
        {
          if (rows[j] == it.row())
          {
            s.soft_constraints_index[j] = index;
          }
        }
      }
      index++;
    }
  }
}

IGL_INLINE void add_soft_constraints_numeric_seperate_final(igl::SLIMData &s)
{
  int v_n = s.v_num;
  int index = 0;
  igl::Timer t;
  t.start();
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      int v_idx = s.b(i);
      s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
    }
  }
  t.stop();
  write_to_file(result_file, "RHS ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  for (int d = 0; d < s.dim; d++)
  {
    for (int i = 0; i < s.b.rows(); i++)
    {
      s.result_vector[s.soft_constraints_index[index]] += s.soft_const_p; // diagonal of matrix
      index++;
    }
  }
  t.stop();
  write_to_file(result_file, "L ADD SOFT CONST", t.getElapsedTimeInMicroSec(), s.first_called);
  t.start();
  tbb::parallel_for(size_t(0), size_t(s.diagonal_index.size()), [&](size_t i) {
    s.result_vector[s.diagonal_index[i]] += s.proximal_p;
  });
  t.stop();
  write_to_file(result_file, "L ADD PROXIMAL", t.getElapsedTimeInMicroSec(), s.first_called);
}

IGL_INLINE double compute_energy(igl::SLIMData &s, Eigen::MatrixXd &V_new)
{
  compute_jacobians(s, V_new);
  return mapping_energy_with_jacobians(s.Ji, s.M, s.slim_energy, s.exp_factor) +
         compute_soft_const_energy(s, s.V, s.F, V_new);
}

IGL_INLINE double compute_soft_const_energy(igl::SLIMData &s,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F,
                                            Eigen::MatrixXd &V_o)
{
  double e = 0;
  for (int i = 0; i < s.b.rows(); i++)
  {
    e += s.soft_const_p * (s.bc.row(i) - V_o.row(s.b(i))).squaredNorm();
  }
  return e;
}

IGL_INLINE void buildRhs(igl::SLIMData &s, const Eigen::SparseMatrix<double> &A)
{
  Eigen::VectorXd f_rhs(s.dim * s.dim * s.f_n);
  f_rhs.setZero();
  if (s.dim == 2)
  {
    /*b = [W11*R11 + W12*R21; (formula (36))
             W11*R12 + W12*R22;
             W21*R11 + W22*R21;
             W21*R12 + W22*R22];*/
    for (int i = 0; i < s.f_n; i++)
    {
      f_rhs(i + 0 * s.f_n) = s.W(i, 0) * s.Ri(i, 0) + s.W(i, 1) * s.Ri(i, 1);
      f_rhs(i + 1 * s.f_n) = s.W(i, 0) * s.Ri(i, 2) + s.W(i, 1) * s.Ri(i, 3);
      f_rhs(i + 2 * s.f_n) = s.W(i, 2) * s.Ri(i, 0) + s.W(i, 3) * s.Ri(i, 1);
      f_rhs(i + 3 * s.f_n) = s.W(i, 2) * s.Ri(i, 2) + s.W(i, 3) * s.Ri(i, 3);
    }
  }
  else
  {
    /*b = [W11*R11 + W12*R21 + W13*R31;
             W11*R12 + W12*R22 + W13*R32;
             W11*R13 + W12*R23 + W13*R33;
             W21*R11 + W22*R21 + W23*R31;
             W21*R12 + W22*R22 + W23*R32;
             W21*R13 + W22*R23 + W23*R33;
             W31*R11 + W32*R21 + W33*R31;
             W31*R12 + W32*R22 + W33*R32;
             W31*R13 + W32*R23 + W33*R33;];*/
    for (int i = 0; i < s.f_n; i++)
    {
      f_rhs(i + 0 * s.f_n) = s.W(i, 0) * s.Ri(i, 0) + s.W(i, 1) * s.Ri(i, 1) + s.W(i, 2) * s.Ri(i, 2);
      f_rhs(i + 1 * s.f_n) = s.W(i, 0) * s.Ri(i, 3) + s.W(i, 1) * s.Ri(i, 4) + s.W(i, 2) * s.Ri(i, 5);
      f_rhs(i + 2 * s.f_n) = s.W(i, 0) * s.Ri(i, 6) + s.W(i, 1) * s.Ri(i, 7) + s.W(i, 2) * s.Ri(i, 8);
      f_rhs(i + 3 * s.f_n) = s.W(i, 3) * s.Ri(i, 0) + s.W(i, 4) * s.Ri(i, 1) + s.W(i, 5) * s.Ri(i, 2);
      f_rhs(i + 4 * s.f_n) = s.W(i, 3) * s.Ri(i, 3) + s.W(i, 4) * s.Ri(i, 4) + s.W(i, 5) * s.Ri(i, 5);
      f_rhs(i + 5 * s.f_n) = s.W(i, 3) * s.Ri(i, 6) + s.W(i, 4) * s.Ri(i, 7) + s.W(i, 5) * s.Ri(i, 8);
      f_rhs(i + 6 * s.f_n) = s.W(i, 6) * s.Ri(i, 0) + s.W(i, 7) * s.Ri(i, 1) + s.W(i, 8) * s.Ri(i, 2);
      f_rhs(i + 7 * s.f_n) = s.W(i, 6) * s.Ri(i, 3) + s.W(i, 7) * s.Ri(i, 4) + s.W(i, 8) * s.Ri(i, 5);
      f_rhs(i + 8 * s.f_n) = s.W(i, 6) * s.Ri(i, 6) + s.W(i, 7) * s.Ri(i, 7) + s.W(i, 8) * s.Ri(i, 8);
    }
  }
  Eigen::VectorXd uv_flat(s.dim * s.v_n);
  for (int i = 0; i < s.dim; i++)
    for (int j = 0; j < s.v_n; j++)
      uv_flat(s.v_n * i + j) = s.V_o(j, i);

  s.rhs = (f_rhs.transpose() * s.WGL_M.asDiagonal() * A).transpose() + s.proximal_p * uv_flat;
}

} // namespace slim
} // namespace igl

IGL_INLINE void igl::slim_update_weights_and_closest_rotations_with_jacobians(const Eigen::MatrixXd &Ji,
                                                                              igl::MappingEnergyType slim_energy,
                                                                              double exp_factor,
                                                                              Eigen::MatrixXd &W,
                                                                              Eigen::MatrixXd &Ri)
{
  const double eps = 1e-8;
  double exp_f = exp_factor;
  const int dim = (Ji.cols() == 4 ? 2 : 3);

  if (dim == 2)
  {
    for (int i = 0; i < Ji.rows(); ++i)
    {
      typedef Eigen::Matrix2d Mat2;
      typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> RMat2;
      typedef Eigen::Vector2d Vec2;
      Mat2 ji, ri, ti, ui, vi;
      Vec2 sing;
      Vec2 closest_sing_vec;
      RMat2 mat_W;
      Vec2 m_sing_new;
      double s1, s2;

      ji(0, 0) = Ji(i, 0);
      ji(0, 1) = Ji(i, 1);
      ji(1, 0) = Ji(i, 2);
      ji(1, 1) = Ji(i, 3);

      igl::polar_svd(ji, ri, ti, ui, sing, vi);

      s1 = sing(0);
      s2 = sing(1);

      // Update Weights according to energy
      switch (slim_energy)
      {
      case igl::MappingEnergyType::ARAP:
      {
        m_sing_new << 1, 1;
        break;
      }
      case igl::MappingEnergyType::SYMMETRIC_DIRICHLET:
      {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case igl::MappingEnergyType::LOG_ARAP:
      {
        double s1_g = 2 * (log(s1) / s1);
        double s2_g = 2 * (log(s2) / s2);
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case igl::MappingEnergyType::CONFORMAL:
      {
        double s1_g = 1 / (2 * s2) - s2 / (2 * pow(s1, 2));
        double s2_g = 1 / (2 * s1) - s1 / (2 * pow(s2, 2));

        double geo_avg = sqrt(s1 * s2);
        double s1_min = geo_avg;
        double s2_min = geo_avg;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
        break;
      }
      case igl::MappingEnergyType::EXP_CONFORMAL:
      {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));

        double geo_avg = sqrt(s1 * s2);
        double s1_min = geo_avg;
        double s2_min = geo_avg;

        double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case igl::MappingEnergyType::EXP_SYMMETRIC_DIRICHLET:
      {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));

        double in_exp = exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      }

      if (std::abs(s1 - 1) < eps)
        m_sing_new(0) = 1;
      if (std::abs(s2 - 1) < eps)
        m_sing_new(1) = 1;
      mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

      W.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 4, Eigen::RowMajor>>(mat_W.data());
      // 2) Update local step (doesn't have to be a rotation, for instance in case of conformal energy)
      Ri.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 4, Eigen::RowMajor>>(ri.data());
    }
  }
  else
  {
    typedef Eigen::Matrix<double, 3, 1> Vec3;
    typedef Eigen::Matrix<double, 3, 3, Eigen::ColMajor> Mat3;
    typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> RMat3;
    Mat3 ji;
    Vec3 m_sing_new;
    Vec3 closest_sing_vec;
    const double sqrt_2 = sqrt(2);
    for (int i = 0; i < Ji.rows(); ++i)
    {
      ji << Ji(i, 0), Ji(i, 1), Ji(i, 2),
          Ji(i, 3), Ji(i, 4), Ji(i, 5),
          Ji(i, 6), Ji(i, 7), Ji(i, 8);

      Mat3 ri, ti, ui, vi;
      Vec3 sing;
      igl::polar_svd(ji, ri, ti, ui, sing, vi);

      double s1 = sing(0);
      double s2 = sing(1);
      double s3 = sing(2);

      // 1) Update Weights
      switch (slim_energy)
      {
      case igl::MappingEnergyType::ARAP:
      {
        m_sing_new << 1, 1, 1;
        break;
      }
      case igl::MappingEnergyType::LOG_ARAP:
      {
        double s1_g = 2 * (log(s1) / s1);
        double s2_g = 2 * (log(s2) / s2);
        double s3_g = 2 * (log(s3) / s3);
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));
        break;
      }
      case igl::MappingEnergyType::SYMMETRIC_DIRICHLET:
      {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        double s3_g = 2 * (s3 - pow(s3, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));
        break;
      }
      case igl::MappingEnergyType::EXP_SYMMETRIC_DIRICHLET:
      {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        double s3_g = 2 * (s3 - pow(s3, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));

        double in_exp = exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;
        s3_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))), sqrt(s3_g / (2 * (s3 - 1)));

        break;
      }
      case igl::MappingEnergyType::CONFORMAL:
      {
        double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

        double s1_g = (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) / common_div;
        double s2_g = (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) / common_div;
        double s3_g = (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) / common_div;

        double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
        double s1_min = closest_s;
        double s2_min = closest_s;
        double s3_min = closest_s;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(s3_g / (2 * (s3 - s3_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min, s3_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
        break;
      }
      case igl::MappingEnergyType::EXP_CONFORMAL:
      {
        // E_conf = (s1^2 + s2^2 + s3^2)/(3*(s1*s2*s3)^(2/3) )
        // dE_conf/ds1 = (-2*(s2*s3)*(s2^2+s3^2 -2*s1^2) ) / (9*(s1*s2*s3)^(5/3))
        // Argmin E_conf(s1): s1 = sqrt(s1^2+s2^2)/sqrt(2)
        double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

        double s1_g = (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) / common_div;
        double s2_g = (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) / common_div;
        double s3_g = (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) / common_div;

        double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow((s1 * s2 * s3), 2. / 3)));
        ;
        double exp_thing = exp(in_exp);

        double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
        double s1_min = closest_s;
        double s2_min = closest_s;
        double s3_min = closest_s;

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;
        s3_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(s3_g / (2 * (s3 - s3_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min, s3_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
      }
      }
      if (std::abs(s1 - 1) < eps)
        m_sing_new(0) = 1;
      if (std::abs(s2 - 1) < eps)
        m_sing_new(1) = 1;
      if (std::abs(s3 - 1) < eps)
        m_sing_new(2) = 1;
      RMat3 mat_W;
      mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

      W.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 9, Eigen::RowMajor>>(mat_W.data());
      // 2) Update closest rotations (not rotations in case of conformal energy)
      Ri.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 9, Eigen::RowMajor>>(ri.data());
    } // for loop end

  } // if dim end
}

IGL_INLINE void igl::slim_buildA(const Eigen::SparseMatrix<double> &Dx,
                                 const Eigen::SparseMatrix<double> &Dy,
                                 const Eigen::SparseMatrix<double> &Dz,
                                 const Eigen::MatrixXd &W,
                                 std::vector<Eigen::Triplet<double>> &IJV)
{
  const int dim = (W.cols() == 4) ? 2 : 3;
  const int f_n = W.rows();
  const int v_n = Dx.cols();

  // formula (35) in paper
  if (dim == 2)
  {
    IJV.reserve(4 * (Dx.outerSize() + Dy.outerSize()));

    /*A = [W11*Dx, W12*Dx;
          W11*Dy, W12*Dy;
          W21*Dx, W22*Dx;
          W21*Dy, W22*Dy];*/
    for (int k = 0; k < Dx.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it)
      {
        int dx_r = it.row();
        int dx_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * W(dx_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(dx_r, v_n + dx_c, val * W(dx_r, 1)));

        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dx_r, dx_c, val * W(dx_r, 2)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dx_r, v_n + dx_c, val * W(dx_r, 3)));
      }
    }

    for (int k = 0; k < Dy.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it)
      {
        int dy_r = it.row();
        int dy_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, dy_c, val * W(dy_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, v_n + dy_c, val * W(dy_r, 1)));

        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dy_r, dy_c, val * W(dy_r, 2)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dy_r, v_n + dy_c, val * W(dy_r, 3)));
      }
    }
  }
  else
  {

    /*A = [W11*Dx, W12*Dx, W13*Dx;
            W11*Dy, W12*Dy, W13*Dy;
            W11*Dz, W12*Dz, W13*Dz;
            W21*Dx, W22*Dx, W23*Dx;
            W21*Dy, W22*Dy, W23*Dy;
            W21*Dz, W22*Dz, W23*Dz;
            W31*Dx, W32*Dx, W33*Dx;
            W31*Dy, W32*Dy, W33*Dy;
            W31*Dz, W32*Dz, W33*Dz;];*/
    IJV.reserve(9 * (Dx.outerSize() + Dy.outerSize() + Dz.outerSize()));
    for (int k = 0; k < Dx.outerSize(); k++)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it)
      {
        int dx_r = it.row();
        int dx_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * W(dx_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(dx_r, v_n + dx_c, val * W(dx_r, 1)));
        IJV.push_back(Eigen::Triplet<double>(dx_r, 2 * v_n + dx_c, val * W(dx_r, 2)));

        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r, dx_c, val * W(dx_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r, v_n + dx_c, val * W(dx_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r, 2 * v_n + dx_c, val * W(dx_r, 5)));

        IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r, dx_c, val * W(dx_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r, v_n + dx_c, val * W(dx_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r, 2 * v_n + dx_c, val * W(dx_r, 8)));
      }
    }

    for (int k = 0; k < Dy.outerSize(); k++)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it)
      {
        int dy_r = it.row();
        int dy_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, dy_c, val * W(dy_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, v_n + dy_c, val * W(dy_r, 1)));
        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, 2 * v_n + dy_c, val * W(dy_r, 2)));

        IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r, dy_c, val * W(dy_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r, v_n + dy_c, val * W(dy_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r, 2 * v_n + dy_c, val * W(dy_r, 5)));

        IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r, dy_c, val * W(dy_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r, v_n + dy_c, val * W(dy_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r, 2 * v_n + dy_c, val * W(dy_r, 8)));
      }
    }

    for (int k = 0; k < Dz.outerSize(); k++)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, k); it; ++it)
      {
        int dz_r = it.row();
        int dz_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r, dz_c, val * W(dz_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r, v_n + dz_c, val * W(dz_r, 1)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r, 2 * v_n + dz_c, val * W(dz_r, 2)));

        IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r, dz_c, val * W(dz_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r, v_n + dz_c, val * W(dz_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r, 2 * v_n + dz_c, val * W(dz_r, 5)));

        IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r, dz_c, val * W(dz_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r, v_n + dz_c, val * W(dz_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r, 2 * v_n + dz_c, val * W(dz_r, 8)));
      }
    }
  }
}
/// Slim Implementation

IGL_INLINE void igl::slim_precompute(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &V_init,
    igl::SLIMData &data,
    igl::MappingEnergyType slim_energy,
    const Eigen::VectorXi &b,
    const Eigen::MatrixXd &bc,
    double soft_p)
{

  data.V = V;
  data.F = F;
  data.V_o = V_init;

  data.v_num = V.rows();
  data.f_num = F.rows();

  data.slim_energy = slim_energy;

  data.b = b;
  data.bc = bc;
  data.soft_const_p = soft_p;

  data.proximal_p = 0.0001;

  igl::doublearea(V, F, data.M);
  data.M /= 2.;
  data.mesh_area = data.M.sum();
  data.mesh_improvement_3d = false; // whether to use a jacobian derived from a real mesh or an abstract regular mesh (used for mesh improvement)
  data.exp_factor = 1.0;            // param used only for exponential energies (e.g exponential symmetric dirichlet)

  assert(F.cols() == 3 || F.cols() == 4);

  igl::slim::pre_calc(data);
  data.energy = igl::slim::compute_energy(data, data.V_o) / data.mesh_area;
}

IGL_INLINE Eigen::MatrixXd igl::slim_solve(igl::SLIMData &data, int iter_num)
{
  for (int i = 0; i < iter_num; i++)
  {
    Eigen::MatrixXd dest_res;
    dest_res = data.V_o;

    // Solve Weighted Proxy
    igl::slim::update_weights_and_closest_rotations(data, dest_res);
    igl::slim::solve_weighted_arap(data, data.V, data.F, dest_res, data.b, data.bc);

    double old_energy = data.energy;

    std::function<double(Eigen::MatrixXd &)> compute_energy = [&](
                                                                  Eigen::MatrixXd &aaa) { return igl::slim::compute_energy(data, aaa); };

    data.energy = igl::flip_avoiding_line_search(data.F, data.V_o, dest_res, compute_energy,
                                                 data.energy * data.mesh_area) /
                  data.mesh_area;
  }
  return data.V_o;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
