#pragma once
#include "igl_inline.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "mkl_symmetric_solver.h"
#include "mkl_routines.hpp"
#include "inline_expansion/NumericExecutor.hpp"

namespace igl
{
struct COTSMOOTHData
{
    // Input
    Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
    Eigen::MatrixXi F; // #F by 3/3 list of mesh faces (triangles/tets)
    double w;
    int method;

    // internal
    Eigen::SparseMatrix<double> Lw;  // weighted L
    Eigen::SparseMatrix<double> L;   // the actual L used in computation
    Eigen::SparseMatrix<double> M;   // mass matrix
    Eigen::SparseMatrix<double> lhs; // left hand side
    Eigen::MatrixXd rhs;             // right hand side

    // output
    Eigen::MatrixXd smoothedV;

    // misc
    bool first_called = true;    // is this the first iteration
    MKLPardisoData pardiso_data; // for the pardiso solver

    // for numeric or mkl
    std::vector<double> result_vector; // store the result, this is essentially the valuePtr
    std::vector<int> L_outer;          // the outerindexptr
    std::vector<int> L_inner;          // the innerindexptr

    // for mkl
    sparse_matrix_t ata_mkl;
    sparse_matrix_t final_result;
    struct matrix_descr symmetric_type;

    // for numeric
    ie::NumericExecutor ex;                 // the holy executor, our lord and savior
    std::vector<std::vector<double>> datas; // datas for executor
    std::vector<int> diagonal_positions;    // to record where the diavonals are in the result_vector
};

IGL_INLINE void constructL(COTSMOOTHData &c);       // construct L from Lw and M
IGL_INLINE void cot_smooth_solve(COTSMOOTHData &c); // the solved vertex will be in smoothedV
IGL_INLINE void build_rhs(COTSMOOTHData &c);        // build the right hand side

} // namespace igl