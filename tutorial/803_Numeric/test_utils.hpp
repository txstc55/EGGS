#pragma once
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tbb/parallel_for.h>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <igl/Timer.h>
#include <map>
#include <set>
#include <fstream>

inline double benchmarkTimer(std::function<void()> op)
{
    igl::Timer t;
    t.start();
    op();
    t.stop();
    return t.getElapsedTimeInMicroSec() / 100.0;
}

#define PROFILE(CODE, MSG)                                                           \
    do                                                                               \
    {                                                                                \
        igl::Timer t;                                                                \
        t.start();                                                                   \
        for (int i = 0; i < 100; i++)                                                \
        {                                                                            \
            (CODE);                                                                  \
        }                                                                            \
        t.stop();                                                                    \
        std::cout << MSG << ": " << t.getElapsedTimeInMicroSec() / 100.0 << " us\n"; \
    } while (0)

// construct sparse matrix from row count, column count, nonzeros, the three arrays
static Eigen::MappedSparseMatrix<double, Eigen::RowMajor> ConstructSparseMatrix(int rowCount, int colCount, int nonZeroCount, double *nonZeroArray, int *rowIndex, int *colIndex)
{
    Eigen::MappedSparseMatrix<double, Eigen::RowMajor> spMap(rowCount, colCount, nonZeroCount, rowIndex, colIndex, nonZeroArray);
    return spMap;
}

static Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ConstructDenseMatrix(int rowCount, int colCount, double *nonZeroArray)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result_matrix(rowCount, colCount);
    unsigned int count = 0;
    for (unsigned int i = 0; i < rowCount; i++)
    {
        for (unsigned int j = 0; j < colCount; j++)
        {
            result_matrix(i, j) = nonZeroArray[count];
            count++;
        }
    }
    return result_matrix;
}

#define PROFILE_EXECUTOR_SINGLE(S, DATA, r_vector, result_file)                                  \
    ;                                                                                            \
    do                                                                                           \
    {                                                                                            \
        auto elapsed_executor_single = benchmarkTimer([&]() {                                    \
            for (int i = 0; i < 100; i++)                                                        \
            {                                                                                    \
                S.ExecuteSingle(DATA, r_vector);                                                 \
            }                                                                                    \
        });                                                                                      \
        std::cout << "NUMERIC EXECUTOR SINGLE THREAD: " << elapsed_executor_single << " us\n";   \
        result_file << "NUMERIC EXECUTOR SINGLE THREAD: " << elapsed_executor_single << " us\n"; \
    } while (0)

#define PROFILE_EXECUTOR_MULTI(S, DATA, r_vector, result_file)                                 \
    ;                                                                                          \
    do                                                                                         \
    {                                                                                          \
        auto elapsed_executor_multi = benchmarkTimer([&]() {                                   \
            for (int i = 0; i < 100; i++)                                                      \
            {                                                                                  \
                S.ExecuteMulti(DATA, r_vector);                                                \
            }                                                                                  \
        });                                                                                    \
        std::cout << "NUMERIC EXECUTOR MULTI THREAD: " << elapsed_executor_multi << " us\n";   \
        result_file << "NUMERIC EXECUTOR MULTI THREAD: " << elapsed_executor_multi << " us\n"; \
    } while (0)

#include <mkl_types.h>
#include <mkl_spblas.h>

template <typename T>
sparse_status_t mkl_create_coo(sparse_matrix_t *A,
                               sparse_index_base_t indexing, MKL_INT num_rows, MKL_INT num_cols,
                               MKL_INT nnz, MKL_INT *rows, MKL_INT *cols, T *values);

/*!
 * Create MKL sparse csr Matrix
 */

template <typename T>
void create_mkl_csr_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor> &M, sparse_matrix_t *A);

template <typename T>
sparse_status_t mkl_export_csr(const sparse_matrix_t A,
                               sparse_index_base_t *indexing, MKL_INT *num_rows, MKL_INT *num_cols,
                               MKL_INT **rows_start, MKL_INT **rows_end, MKL_INT **col_indx,
                               T **values);

#define PROFILE_MKL_SINGLE_PIPE(MATRIX_VECTOR, R, result_file)                                                                                                                                                                          \
    do                                                                                                                                                                                                                                  \
    {                                                                                                                                                                                                                                   \
        sparse_status_t status;                                                                                                                                                                                                         \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                                                                                    \
        struct matrix_descr general_type;                                                                                                                                                                                               \
        general_type.type = SPARSE_MATRIX_TYPE_GENERAL;                                                                                                                                                                                 \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                         \
        {                                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                                        \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                                                                                  \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                                                                              \
            matrix_vectors.push_back(tmp);                                                                                                                                                                                              \
        }                                                                                                                                                                                                                               \
        for (int i = 0; i < MATRIX_VECTOR.size() - 1; i++)                                                                                                                                                                              \
        {                                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                                        \
            result_vectors.push_back(tmp);                                                                                                                                                                                              \
        }                                                                                                                                                                                                                               \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_NNZ_COUNT, &result_vectors[0]);                         \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[0]);              \
        for (unsigned int i = 2; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                         \
        {                                                                                                                                                                                                                               \
            status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_NNZ_COUNT, &result_vectors[i - 1]);             \
            status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[i - 1]);  \
        }                                                                                                                                                                                                                               \
        mkl_set_num_threads_local(1);                                                                                                                                                                                                   \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                                                                                           \
            for (int j = 0; j < 100; j++)                                                                                                                                                                                               \
            {                                                                                                                                                                                                                           \
                status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[0]);             \
                for (unsigned int i = 2; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                 \
                {                                                                                                                                                                                                                       \
                    status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[i - 1]); \
                }                                                                                                                                                                                                                       \
            }                                                                                                                                                                                                                           \
        });                                                                                                                                                                                                                             \
        std::cout << "MKL_SINGLE " << MATRIX_VECTOR.size() << " MATRICES: " << elapsed << " us\n";                                                                                                                                      \
        result_file << "MKL_SINGLE " << MATRIX_VECTOR.size() << " MATRICES: " << elapsed << " us\n";                                                                                                                                    \
        R = result_vectors[result_vectors.size() - 1];                                                                                                                                                                                  \
        mkl_free_buffers();                                                                                                                                                                                                             \
    } while (0)

#define PROFILE_MKL_MULTI_PIPE(MATRIX_VECTOR, R, result_file)                                                                                                                                                                           \
    do                                                                                                                                                                                                                                  \
    {                                                                                                                                                                                                                                   \
        sparse_status_t status;                                                                                                                                                                                                         \
        struct matrix_descr general_type;                                                                                                                                                                                               \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                                                                                    \
        general_type.type = SPARSE_MATRIX_TYPE_GENERAL;                                                                                                                                                                                 \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                         \
        {                                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                                        \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                                                                                  \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                                                                              \
            matrix_vectors.push_back(tmp);                                                                                                                                                                                              \
        }                                                                                                                                                                                                                               \
        for (int i = 0; i < MATRIX_VECTOR.size() - 1; i++)                                                                                                                                                                              \
        {                                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                                        \
            result_vectors.push_back(tmp);                                                                                                                                                                                              \
        }                                                                                                                                                                                                                               \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_NNZ_COUNT, &result_vectors[0]);                         \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[0]);              \
        for (unsigned int i = 2; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                         \
        {                                                                                                                                                                                                                               \
            status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_NNZ_COUNT, &result_vectors[i - 1]);             \
            status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[i - 1]);  \
        }                                                                                                                                                                                                                               \
        mkl_set_num_threads_local(0);                                                                                                                                                                                                   \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                                                                                           \
            for (int j = 0; j < 100; j++)                                                                                                                                                                                               \
            {                                                                                                                                                                                                                           \
                status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[0]);             \
                for (unsigned int i = 2; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                                 \
                {                                                                                                                                                                                                                       \
                    status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[i - 2], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[i], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[i - 1]); \
                }                                                                                                                                                                                                                       \
            }                                                                                                                                                                                                                           \
        });                                                                                                                                                                                                                             \
        std::cout << "MKL_MULTI " << MATRIX_VECTOR.size() << " MATRICES: " << elapsed << " us\n";                                                                                                                                       \
        result_file << "MKL_MULTI " << MATRIX_VECTOR.size() << " MATRICES: " << elapsed << " us\n";                                                                                                                                     \
        R = result_vectors[result_vectors.size() - 1];                                                                                                                                                                                  \
        mkl_free_buffers();                                                                                                                                                                                                             \
    } while (0)

#define PROFILE_MKL_CONST_OP_SINGLE(MATRIX_VECTOR, R, result_file)                                                                                                                                                      \
    do                                                                                                                                                                                                                  \
    {                                                                                                                                                                                                                   \
        sparse_status_t status;                                                                                                                                                                                         \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                                                                    \
        result_vectors.resize(3);                                                                                                                                                                                       \
        struct matrix_descr general_type;                                                                                                                                                                               \
        general_type.type = SPARSE_MATRIX_TYPE_GENERAL;                                                                                                                                                                 \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                         \
        {                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                        \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                                                                  \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                                                              \
            matrix_vectors.push_back(tmp);                                                                                                                                                                              \
        }                                                                                                                                                                                                               \
        mkl_set_num_threads_local(1);                                                                                                                                                                                   \
        auto elapsed_add = benchmarkTimer([&]() {                                                                                                                                                                       \
            for (int j = 0; j < 100; j++)                                                                                                                                                                               \
            {                                                                                                                                                                                                           \
                status = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, matrix_vectors[0], 3.78, matrix_vectors[1], &result_vectors[0]);                                                                              \
                status = mkl_sparse_d_add(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[1], 6.942, matrix_vectors[2], &result_vectors[1]);                                                                                 \
            }                                                                                                                                                                                                           \
        });                                                                                                                                                                                                             \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_NNZ_COUNT, &result_vectors[2]);             \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[2]);  \
        auto elapsed_mult = benchmarkTimer([&]() {                                                                                                                                                                      \
            for (int j = 0; j < 100; j++)                                                                                                                                                                               \
            {                                                                                                                                                                                                           \
                status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[2]); \
            }                                                                                                                                                                                                           \
        });                                                                                                                                                                                                             \
        std::cout << "MKL CONST TEST SINGLE: " << elapsed_add + elapsed_mult << " us\n";                                                                                                                                \
        result_file << "MKL CONST TEST SINGLE: " << elapsed_add + elapsed_mult << " us\n";                                                                                                                              \
        R = result_vectors[2];                                                                                                                                                                                          \
        mkl_free_buffers();                                                                                                                                                                                             \
    } while (0)

#define PROFILE_MKL_CONST_OP_MULT(MATRIX_VECTOR, R, result_file)                                                                                                                                                        \
    do                                                                                                                                                                                                                  \
    {                                                                                                                                                                                                                   \
        sparse_status_t status;                                                                                                                                                                                         \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                                                                    \
        struct matrix_descr general_type;                                                                                                                                                                               \
        for (int i = 0; i < 3; i++)                                                                                                                                                                                     \
        {                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                        \
            result_vectors.push_back(tmp);                                                                                                                                                                              \
        }                                                                                                                                                                                                               \
        general_type.type = SPARSE_MATRIX_TYPE_GENERAL;                                                                                                                                                                 \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                                                                         \
        {                                                                                                                                                                                                               \
            sparse_matrix_t tmp;                                                                                                                                                                                        \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                                                                  \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                                                              \
            matrix_vectors.push_back(tmp);                                                                                                                                                                              \
        }                                                                                                                                                                                                               \
        mkl_set_num_threads_local(0);                                                                                                                                                                                   \
        auto elapsed_add = benchmarkTimer([&]() {                                                                                                                                                                       \
            for (int j = 0; j < 100; j++)                                                                                                                                                                               \
            {                                                                                                                                                                                                           \
                status = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, matrix_vectors[0], 3.78, matrix_vectors[1], &result_vectors[0]);                                                                              \
                status = mkl_sparse_d_add(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[1], 6.942, matrix_vectors[2], &result_vectors[1]);                                                                                 \
            }                                                                                                                                                                                                           \
        });                                                                                                                                                                                                             \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_NNZ_COUNT, &result_vectors[2]);             \
        status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &result_vectors[2]);  \
        auto elapsed_mult = benchmarkTimer([&]() {                                                                                                                                                                      \
            for (int j = 0; j < 100; j++)                                                                                                                                                                               \
            {                                                                                                                                                                                                           \
                status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, general_type, result_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, result_vectors[1], SPARSE_STAGE_FINALIZE_MULT, &result_vectors[2]); \
            }                                                                                                                                                                                                           \
        });                                                                                                                                                                                                             \
        std::cout << "MKL CONST TEST MULT: " << elapsed_add + elapsed_mult << " us\n";                                                                                                                                  \
        result_file << "MKL CONST TEST MULT: " << elapsed_add + elapsed_mult << " us\n";                                                                                                                                \
        R = result_vectors[2];                                                                                                                                                                                          \
        mkl_free_buffers();                                                                                                                                                                                             \
    } while (0)

#define PROFILE_MKL_SINGLE_SYPR(MATRIX_VECTOR, R, result_file)                                                                                                              \
    do                                                                                                                                                                      \
    {                                                                                                                                                                       \
        sparse_status_t status;                                                                                                                                             \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                        \
        struct matrix_descr symmetric_type;                                                                                                                                 \
        symmetric_type.type = SPARSE_MATRIX_TYPE_SYMMETRIC;                                                                                                                 \
        symmetric_type.mode = SPARSE_FILL_MODE_UPPER;                                                                                                                       \
        symmetric_type.diag = SPARSE_DIAG_NON_UNIT;                                                                                                                         \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                             \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                      \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                  \
            matrix_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        for (int i = 0; i < MATRIX_VECTOR.size() - 1; i++)                                                                                                                  \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            result_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_NNZ_COUNT);             \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT_NO_VAL);  \
        mkl_set_num_threads_local(1);                                                                                                                                       \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                               \
            for (int j = 0; j < 100; j++)                                                                                                                                   \
            {                                                                                                                                                               \
                status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT); \
            }                                                                                                                                                               \
        });                                                                                                                                                                 \
        std::cout << "MKL_SINGLE SYPR: " << elapsed << " us\n";                                                                                                             \
        result_file << "MKL_SINGLE SYPR: " << elapsed << " us\n";                                                                                                           \
        R = result_vectors[0];                                                                                                                                              \
        mkl_free_buffers();                                                                                                                                                 \
    } while (0)

#define PROFILE_MKL_MULTI_SYPR(MATRIX_VECTOR, R, result_file)                                                                                                               \
    do                                                                                                                                                                      \
    {                                                                                                                                                                       \
        sparse_status_t status;                                                                                                                                             \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                        \
        struct matrix_descr symmetric_type;                                                                                                                                 \
        symmetric_type.type = SPARSE_MATRIX_TYPE_SYMMETRIC;                                                                                                                 \
        symmetric_type.mode = SPARSE_FILL_MODE_UPPER;                                                                                                                       \
        symmetric_type.diag = SPARSE_DIAG_NON_UNIT;                                                                                                                         \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                             \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                      \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                  \
            matrix_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        for (int i = 0; i < MATRIX_VECTOR.size() - 1; i++)                                                                                                                  \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            result_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_NNZ_COUNT);             \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT_NO_VAL);  \
        mkl_set_num_threads_local(0);                                                                                                                                       \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                               \
            for (int j = 0; j < 100; j++)                                                                                                                                   \
            {                                                                                                                                                               \
                status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT); \
            }                                                                                                                                                               \
        });                                                                                                                                                                 \
        std::cout << "MKL_MULTI SYPR: " << elapsed << " us\n";                                                                                                              \
        result_file << "MKL_SINGLE SYPR: " << elapsed << " us\n";                                                                                                           \
        R = result_vectors[0];                                                                                                                                              \
        mkl_free_buffers();                                                                                                                                                 \
    } while (0)

#define PROFILE_MKL_SINGLE_SYRK(MATRIX_VECTOR, R, result_file)                                                                                                              \
    do                                                                                                                                                                      \
    {                                                                                                                                                                       \
        sparse_status_t status;                                                                                                                                             \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                        \
        struct matrix_descr symmetric_type;                                                                                                                                 \
        symmetric_type.type = SPARSE_MATRIX_TYPE_DIAGONAL;                                                                                                                  \
        symmetric_type.diag = SPARSE_DIAG_UNIT;                                                                                                                             \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                             \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                      \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                  \
            matrix_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        sparse_matrix_t tmp;                                                                                                                                                \
        matrix_vectors.push_back(tmp);                                                                                                                                      \
        for (int i = 0; i < 1; i++)                                                                                                                                         \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            result_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_NNZ_COUNT);             \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT_NO_VAL);  \
        mkl_set_num_threads_local(1);                                                                                                                                       \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                               \
            for (int j = 0; j < 100; j++)                                                                                                                                   \
            {                                                                                                                                                               \
                status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT); \
            }                                                                                                                                                               \
        });                                                                                                                                                                 \
        std::cout << "MKL_SINGLE SYRK: " << elapsed << " us\n";                                                                                                             \
        result_file << "MKL_SINGLE SYPR: " << elapsed << " us\n";                                                                                                           \
        R = result_vectors[0];                                                                                                                                              \
        mkl_free_buffers();                                                                                                                                                 \
    } while (0)

#define PROFILE_MKL_MULTI_SYRK(MATRIX_VECTOR, R, result_file)                                                                                                               \
    do                                                                                                                                                                      \
    {                                                                                                                                                                       \
        sparse_status_t status;                                                                                                                                             \
        std::vector<sparse_matrix_t> matrix_vectors, result_vectors;                                                                                                        \
        struct matrix_descr symmetric_type;                                                                                                                                 \
        symmetric_type.type = SPARSE_MATRIX_TYPE_DIAGONAL;                                                                                                                  \
        symmetric_type.diag = SPARSE_DIAG_UNIT;                                                                                                                             \
        for (unsigned int i = 0; i < MATRIX_VECTOR.size(); i++)                                                                                                             \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            MKL_INT nnz = MATRIX_VECTOR[i].nonZeros();                                                                                                                      \
            create_mkl_csr_matrix(MATRIX_VECTOR[i], &tmp);                                                                                                                  \
            matrix_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        sparse_matrix_t tmp;                                                                                                                                                \
        matrix_vectors.push_back(tmp);                                                                                                                                      \
        for (int i = 0; i < 1; i++)                                                                                                                                         \
        {                                                                                                                                                                   \
            sparse_matrix_t tmp;                                                                                                                                            \
            result_vectors.push_back(tmp);                                                                                                                                  \
        }                                                                                                                                                                   \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_NNZ_COUNT);             \
        status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT_NO_VAL);  \
        mkl_set_num_threads_local(0);                                                                                                                                       \
        auto elapsed = benchmarkTimer([&]() {                                                                                                                               \
            for (int j = 0; j < 100; j++)                                                                                                                                   \
            {                                                                                                                                                               \
                status = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE, matrix_vectors[0], matrix_vectors[1], symmetric_type, &result_vectors[0], SPARSE_STAGE_FINALIZE_MULT); \
            }                                                                                                                                                               \
        });                                                                                                                                                                 \
        std::cout << "MKL_MULTI SYRK: " << elapsed << " us\n";                                                                                                              \
        result_file << "MKL_MULTI SYRK: " << elapsed << " us\n";                                                                                                            \
        R = result_vectors[0];                                                                                                                                              \
        mkl_free_buffers();                                                                                                                                                 \
    } while (0)

#define PROFILE_MKL_MULTI_SPMV(MATRIX, x, y, result_file)                                              \
    do                                                                                                 \
    {                                                                                                  \
        sparse_status_t status;                                                                        \
        struct matrix_descr general_type;                                                              \
        sparse_matrix_t A;                                                                             \
        general_type.type = SPARSE_MATRIX_TYPE_GENERAL;                                                \
        create_mkl_csr_matrix(MATRIX, &A);                                                             \
        mkl_set_num_threads_local(0);                                                                  \
        auto elapsed = benchmarkTimer([&]() {                                                          \
            for (int j = 0; j < 100; j++)                                                              \
            {                                                                                          \
                status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, general_type, x, 0, y); \
            }                                                                                          \
        });                                                                                            \
        std::cout << "MKL_MULTI SPMV: " << elapsed << " us\n";                                         \
        result_file << "MKL_MULTI SPMV: " << elapsed << " us\n";                                       \
        mkl_free_buffers();                                                                            \
    } while (0)

template <typename T>
static std::tuple<int, int, int, std::vector<int>, std::vector<int>, std::vector<T>> extract_value(sparse_matrix_t &mkl_matrix)
{
    sparse_status_t status;
    sparse_index_base_t indexing;
    MKL_INT num_rows, num_cols;
    MKL_INT *row_start, *row_end, *col_indx;
    T *values;
    status = mkl_export_csr<T>(mkl_matrix, &indexing, &num_rows, &num_cols,
                               &row_start, &row_end,
                               &col_indx, &values);
    int count = 0;
    std::vector<int> outerindex;
    std::vector<int> innerindex;
    std::vector<T> valueptr;
    outerindex.reserve(num_rows + 1);
    innerindex.reserve(num_rows);
    valueptr.reserve(num_rows);
    for (int i = 0; i < num_rows; ++i)
    {
        outerindex.push_back(row_start[i]);
        std::vector<std::pair<int, T>> sorted_index;
        for (int j = row_start[i]; j < row_end[i]; ++j)
        {
            sorted_index.push_back({col_indx[j], values[count]});
            ++count;
        }
        sort(sorted_index.begin(), sorted_index.end());
        for (auto p : sorted_index)
        {
            innerindex.push_back(std::get<0>(p));
            valueptr.push_back(std::get<1>(p));
        }
    }
    outerindex.push_back(count);
    return {num_rows, num_cols, count, outerindex, innerindex, valueptr};
}

static Eigen::SparseMatrix<double, Eigen::RowMajor> read_mtx(std::string filename)
{
    std::ifstream stream(filename);
    int row, col, nnz;
    double val;
    std::string formats;
    getline(stream, formats);
    bool if_symmetric = formats.find("symmetric") != std::string::npos;

    std::string comments;
    while (true)
    {
        getline(stream, comments);
        if (comments[0] != '%')
        {
            std::stringstream specs(comments);
            specs >> row >> col >> nnz;
            if (if_symmetric)
            {
                nnz = nnz * 2;
            }
            break;
        }
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(row, col);
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(nnz);
    while (stream >> row >> col >> val)
    {
        trip.push_back(Eigen::Triplet<double>(row - 1, col - 1, val));
        if (if_symmetric)
        {
            trip.push_back(Eigen::Triplet<double>(col - 1, row - 1, val));
        }
    }
    M.setFromTriplets(trip.begin(), trip.end());
    M.makeCompressed();
    return M;
}

static void prettify_test_name(std::string info)
{
    std::cout << "========================================================================\n";
    std::string information = "|| " + info;
    information += std::string(70 - information.size(), ' ') + "||\n";
    std::cout << information << "========================================================================\n";
}
