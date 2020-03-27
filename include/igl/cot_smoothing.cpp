#include "cot_smoothing.h"
#include "cotmatrix.h"            // for getting the cotangent matrix Lw
#include "massmatrix.h"           // for getting the mass matrix M
#include "mkl_symmetric_solver.h" // for solving the matrix
#include <iostream>
#include "Timer.h"
#include "test_record_util.h"

#include "inline_expansion/utils.hpp"
#include "ie_helper.hpp"
#include <tbb/parallel_for.h>
#include "mkl.h"
using namespace Eigen;
using namespace std;
namespace igl
{

IGL_INLINE void build_linear_system_eigen(COTSMOOTHData &c);            // build by eigen doing AtBA + cD
IGL_INLINE void build_linear_system_numeric_together(COTSMOOTHData &c); // build by inline-expansion
IGL_INLINE void build_linear_system_mkl(COTSMOOTHData &c);              // build by mkl
IGL_INLINE void build_linear_system_numeric_seperate(COTSMOOTHData &c); // build by inline-expansion, but do ata and the plus terms seperate

std::ofstream result_file;

IGL_INLINE void constructL(COTSMOOTHData &c)
{
    // determine if this is the first time of calling
    if (c.first_called)
    {
        c.smoothedV = c.V;         // copy the vertex and never touch V again
        cotmatrix(c.V, c.F, c.Lw); // compute the cotmatrix once
        c.Lw.makeCompressed();     // compress this
        c.L = c.Lw;                // just allocating the nonzeros
    }
    massmatrix(c.smoothedV, c.F, MASSMATRIX_TYPE_VORONOI, c.M); // getting the mass matrix, which is diagonal
    c.M.makeCompressed();                                       // compress it once

    // we need to construct L each time the mass matrix is changed
    // since L = M.inverse()*Lw
    int ind = 0;
    for (int i = 0; i < c.L.outerSize(); i++)
    {
        for (SparseMatrix<double>::InnerIterator it(c.L, i); it; ++it)
        {
            c.L.valuePtr()[ind] = c.Lw.valuePtr()[ind] / c.M.valuePtr()[it.row()];
            ind++;
        }
    }
    c.L.makeCompressed();
    // by this time we have finished construction of L and M
}

IGL_INLINE void build_rhs(COTSMOOTHData &c)
{
    // if first time init
    if (c.first_called)
    {
        c.rhs.resize(c.V.rows(), 3);
    }

    // this is computing rhs = wMp
    for (int i = 0; i < c.smoothedV.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            c.rhs(i, j) = c.smoothedV(i, j) * c.M.valuePtr()[i] * c.w;
        }
    }
}

IGL_INLINE void build_linear_system_eigen(COTSMOOTHData &c)
{
    result_file << "START BUILDING LINEAR SYSTEM USING EIGEN\n";
    igl::Timer t;
    t.start();
    constructL(c); // construct L
    t.stop();
    write_to_file(result_file, "CONSTRUCTING L", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    build_rhs(c);
    t.stop();
    write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    Eigen::SparseMatrix<double> Lt = c.L.transpose();
    t.stop();
    write_to_file(result_file, "TAKING TRANSPOSE", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    c.lhs = Lt * c.M * c.L + c.w * c.M;
    t.stop();
    write_to_file(result_file, "EIGEN COMPUTATION", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    c.lhs = c.lhs.triangularView<Eigen::Lower>(); // because pardiso only wants the triangular
    c.lhs.makeCompressed();
    t.stop();
    write_to_file(result_file, "GETTING TRIANGULAR", t.getElapsedTimeInMicroSec(), c.first_called);
}

IGL_INLINE void build_linear_system_mkl(COTSMOOTHData &c)
{
    result_file << "START BUILDING LINEAR SYSTEM USING MKL\n";
    igl::Timer t;
    t.start();
    constructL(c); // construct L
    t.stop();
    write_to_file(result_file, "CONSTRUCTING L", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    build_rhs(c);
    t.stop();
    write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    sparse_matrix_t L_mkl;
    create_mkl_csr_matrix(c.L, &L_mkl);
    t.stop();
    write_to_file(result_file, "CONVERTING L TO MKL FORMAT", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    sparse_matrix_t M_mkl;
    create_mkl_csr_matrix(c.M, &M_mkl);
    t.stop();
    write_to_file(result_file, "CONVERTING M TO MKL FORMAT", t.getElapsedTimeInMicroSec(), c.first_called);
    if (c.first_called)
    {
        mkl_sypr_pre(L_mkl, M_mkl, c.symmetric_type, &c.ata_mkl);
    }
    t.start();
    mkl_sypr_final(L_mkl, M_mkl, c.symmetric_type, &c.ata_mkl);
    t.stop();
    write_to_file(result_file, "MKL ATA", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    mkl_sparse_add(M_mkl, c.ata_mkl, c.w, &c.final_result);
    t.stop();
    write_to_file(result_file, "MKL ADDING WEIGHTED M", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    if (c.first_called)
    {
        export_csr_from_mkl(c.final_result, c.L_outer, c.L_inner, c.result_vector);
    }
    else
    {
        export_val_from_mkl(c.final_result, c.result_vector);
    }
    t.stop();
    write_to_file(result_file, "EXPORTING TO CSR", t.getElapsedTimeInMicroSec(), c.first_called);
}

IGL_INLINE void build_linear_system_numeric_together(COTSMOOTHData &c)
{
    result_file << "START BUILDING LINEAR SYSTEM USING NUMERIC TOGETHER\n";
    igl::Timer t;
    t.start();
    constructL(c); // construct L
    t.stop();
    write_to_file(result_file, "CONSTRUCTING L", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    build_rhs(c);
    t.stop();
    write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), c.first_called);
    if (c.first_called)
    {
        c.datas.resize(2);
        c.datas[0].resize(c.L.nonZeros());
        c.datas[1].resize(c.L.rows());
        SparseMatrix<ie::NumericType> L_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(c.L, 0);
        SparseMatrix<ie::NumericType> M_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(c.M, 1);
        SparseMatrix<ie::NumericType> result_numeric = Eigen::SparseMatrix<ie::NumericType>(L_numeric.transpose()) * M_numeric * L_numeric + c.w * M_numeric;
        result_numeric = result_numeric.triangularView<Eigen::Lower>();
        c.ex = ie::NumericExecutor(result_numeric, 0);
        c.result_vector.resize(result_numeric.nonZeros(), 0);
        // copy the outer index pointer and inner index pointer to s
        c.L_outer.resize(result_numeric.rows() + 1);
        c.L_inner.resize(result_numeric.nonZeros());
        c.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
        c.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
    }
    t.start();
    c.datas[0].assign(c.L.valuePtr(), c.L.valuePtr() + c.L.nonZeros());
    c.datas[1].assign(c.M.valuePtr(), c.M.valuePtr() + c.M.nonZeros());
    t.stop();
    write_to_file(result_file, "ASSIGNING DATAS", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    c.ex.ExecuteMulti(c.datas, c.result_vector);
    t.stop();
    write_to_file(result_file, "COMPUTE EVERYTHING", t.getElapsedTimeInMicroSec(), c.first_called);
}

IGL_INLINE void build_linear_system_numeric_seperate(COTSMOOTHData &c)
{
    result_file << "START BUILDING LINEAR SYSTEM USING NUMERIC SEPERATE\n";
    igl::Timer t;
    t.start();
    constructL(c); // construct L
    t.stop();
    write_to_file(result_file, "CONSTRUCTING L", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    build_rhs(c);
    t.stop();
    write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), c.first_called);
    if (c.first_called)
    {
        c.datas.resize(2);
        c.datas[0].resize(c.L.nonZeros());
        c.datas[1].resize(c.L.rows());
        SparseMatrix<ie::NumericType> L_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(c.L, 0);
        SparseMatrix<ie::NumericType> M_numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(c.M, 1);
        SparseMatrix<ie::NumericType> result_numeric = Eigen::SparseMatrix<ie::NumericType>(L_numeric.transpose()) * M_numeric * L_numeric;
        result_numeric = result_numeric.triangularView<Eigen::Lower>();
        c.ex = ie::NumericExecutor(result_numeric, 0);
        c.result_vector.resize(result_numeric.nonZeros(), 0);
        // copy the outer index pointer and inner index pointer to s
        c.L_outer.resize(result_numeric.rows() + 1);
        c.L_inner.resize(result_numeric.nonZeros());
        c.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
        c.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
        int index = 0;
        c.diagonal_positions.reserve(c.L.rows());
        for (int i = 0; i < result_numeric.outerSize(); i++)
        {
            for (SparseMatrix<ie::NumericType>::InnerIterator it(result_numeric, i); it; ++it)
            {
                if (it.row() == it.col())
                {
                    c.diagonal_positions.push_back(index);
                }
                index++;
            }
        }
    }
    t.start();
    c.datas[0].assign(c.L.valuePtr(), c.L.valuePtr() + c.L.nonZeros());
    c.datas[1].assign(c.M.valuePtr(), c.M.valuePtr() + c.M.nonZeros());
    t.stop();
    write_to_file(result_file, "ASSIGNING DATAS", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    c.ex.ExecuteMulti(c.datas, c.result_vector);
    t.stop();
    write_to_file(result_file, "COMPUTE ATDA", t.getElapsedTimeInMicroSec(), c.first_called);
    t.start();
    tbb::parallel_for(size_t(0), size_t(c.L.rows()), [&](size_t i) {
        c.result_vector[c.diagonal_positions[i]] += c.w * c.datas[1][i];
    });
    t.stop();
    write_to_file(result_file, "ADDING WEIGHTED M", t.getElapsedTimeInMicroSec(), c.first_called);
}

IGL_INLINE void cot_smooth_solve(COTSMOOTHData &c)
{
    igl::Timer t;
    if (c.first_called)
    {
        result_file.open("result_cot.txt");
    }
    switch (c.method)
    {
    case 0:
        build_linear_system_eigen(c);
        break;
    case 2:
        build_linear_system_mkl(c);
        break;
    case 3:
        build_linear_system_numeric_together(c);
        break;
    case 4:
        build_linear_system_numeric_seperate(c);
        break;
    default:
        build_linear_system_eigen(c);
        break;
    }
    std::vector<double> solved(c.V.rows() * 3);
    // support the lhs, if first called, do symbolic factorization
    if (c.first_called)
    {
        pardiso_init(c.pardiso_data);
        if (c.result_vector.size() == 0)
        {
            pardiso_support_matrix(c.pardiso_data, c.lhs);
        }
        else
        {
            pardiso_support_matrix(c.pardiso_data, c.L_outer.data(), c.L_inner.data(), c.result_vector.data(), c.V.rows());
        }
        pardiso_symbolic_factor(c.pardiso_data);
    }
    else
    {
        if (c.result_vector.size() == 0)
        {
            pardiso_support_value(c.pardiso_data, c.lhs.valuePtr());
        }
        else
        {
            pardiso_support_value(c.pardiso_data, c.result_vector.data());
        }
    }
    t.start();
    pardiso_numeric_factor(c.pardiso_data);
    t.stop();
    write_to_file(result_file, "MKL NUMERIC FACTORIZATION", t.getElapsedTimeInMicroSec(), c.first_called);
    // solve for x
    t.start();
    pardiso_solve(c.pardiso_data, solved.data(), c.rhs.data());
    // solve for y
    pardiso_solve(c.pardiso_data, solved.data() + c.rhs.rows(), c.rhs.data() + c.rhs.rows());
    // solve for z
    pardiso_solve(c.pardiso_data, solved.data() + c.rhs.rows() * 2, c.rhs.data() + 2 * c.rhs.rows());
    t.stop();
    write_to_file(result_file, "MKL SOLVING", t.getElapsedTimeInMicroSec(), c.first_called);

    // putting it back into the smoothedV
    c.smoothedV = Map<MatrixXd>(solved.data(), c.rhs.rows(), 3);
    c.first_called = false;
}

} // namespace igl