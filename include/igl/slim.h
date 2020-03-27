// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SLIM_H
#define SLIM_H

#include "igl_inline.h"
#include "MappingEnergyType.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mkl_symmetric_solver.h"
#include "mkl_routines.hpp"
#include "inline_expansion/NumericExecutor.hpp"

// This option makes the iterations faster (all except the first) by caching the
// sparsity pattern of the matrix involved in the assembly. It should be on if you plan to do many iterations, off if you have to change the matrix structure at every iteration.
#define SLIM_CACHED

#ifdef SLIM_CACHED
#include <igl/AtA_cached.h>
#endif

namespace igl
{

// Compute a SLIM map as derived in "Scalable Locally Injective Maps" [Rabinovich et al. 2016].
struct SLIMData
{
  // Input
  Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
  Eigen::MatrixXi F; // #F by 3/3 list of mesh faces (triangles/tets)
  MappingEnergyType slim_energy;

  // Optional Input
  // soft constraints
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  double soft_const_p;

  double exp_factor;        // used for exponential energies, ignored otherwise
  bool mesh_improvement_3d; // only supported for 3d

  // Output
  Eigen::MatrixXd V_o; // #V by dim list of mesh vertex positions (dim = 2 for parametrization, 3 otherwise)
  double energy;       // objective value

  // INTERNAL
  Eigen::VectorXd M;
  double mesh_area;
  double avg_edge_length;
  int v_num;
  int f_num;
  double proximal_p;

  Eigen::VectorXd WGL_M;
  Eigen::VectorXd rhs;
  Eigen::MatrixXd Ri, Ji;
  Eigen::MatrixXd W;
  Eigen::SparseMatrix<double> Dx, Dy, Dz;
  int f_n, v_n;
  bool first_solve;
  bool has_pre_calc = false;
  int dim;

#ifdef SLIM_CACHED
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXi A_data;
  Eigen::SparseMatrix<double> AtA;
  igl::AtA_cached_data AtA_data;
#endif

  // using mkl solver
  MKLPardisoData pardiso_data;
  bool first_called = true;

  // using inline-expansion
  ie::NumericExecutor ex = ie::NumericExecutor();               // the executor itself
  std::vector<double> result_vector;                            // the vector that holds the values when the calculation is done
  std::vector<int> soft_constraints_pos;                        // record how we order the soft constraints, it is basically a map except the key is the position itself
  std::vector<std::vector<double>> datas;                       // the values used for executor
  std::vector<Eigen::Triplet<double>> soft_constraints_triplet; // will only be used once, construct the soft constraints matrix
  int n_rows;                                                   // the number of rows of the AtbA
  std::vector<int> L_outer;                                     // the outerindexptr
  std::vector<int> L_inner;                                     // the innerindexptr
  std::vector<int> A_triplet_order;                             // record how we insert the triplet's value, basically a map except the key is the position itself

  // using inline-expansion, but seperate the steps
  std::vector<int> soft_constraints_index; // record where at the result_vector we need to add the soft constraints
  std::vector<int> diagonal_index;         // record where the diagonal entries are in the result_vector

  // using mkl
  sparse_matrix_t ata_mkl;
  sparse_matrix_t final_result;
  struct matrix_descr symmetric_type;



  // specify the method type
  int method_type;
};

// Compute necessary information to start using SLIM
// Inputs:
//      V           #V by 3 list of mesh vertex positions
//      F           #F by 3/3 list of mesh faces (triangles/tets)
//    b           list of boundary indices into V
//    bc          #b by dim list of boundary conditions
//    soft_p      Soft penalty factor (can be zero)
//    slim_energy Energy to minimize
IGL_INLINE void slim_precompute(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &V_init,
    SLIMData &data,
    MappingEnergyType slim_energy,
    const Eigen::VectorXi &b,
    const Eigen::MatrixXd &bc,
    double soft_p);

// Run iter_num iterations of SLIM
// Outputs:
//    V_o (in SLIMData): #V by dim list of mesh vertex positions
IGL_INLINE Eigen::MatrixXd slim_solve(SLIMData &data, int iter_num);

// Internal Routine. Exposed for Integration with SCAF
IGL_INLINE void slim_update_weights_and_closest_rotations_with_jacobians(const Eigen::MatrixXd &Ji,
                                                                         igl::MappingEnergyType slim_energy,
                                                                         double exp_factor,
                                                                         Eigen::MatrixXd &W,
                                                                         Eigen::MatrixXd &Ri);

IGL_INLINE void slim_buildA(const Eigen::SparseMatrix<double> &Dx,
                            const Eigen::SparseMatrix<double> &Dy,
                            const Eigen::SparseMatrix<double> &Dz,
                            const Eigen::MatrixXd &W,
                            std::vector<Eigen::Triplet<double>> &IJV);
} // namespace igl

#ifndef IGL_STATIC_LIBRARY
#include "slim.cpp"
#endif

#endif // SLIM_H
