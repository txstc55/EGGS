#include "mkl_routines.hpp"
#include "mkl.h"
#include <iostream>

namespace igl
{
IGL_INLINE void create_mkl_csr_matrix(Eigen::SparseMatrix<double, Eigen::ColMajor> &M, sparse_matrix_t *A)
{
    sparse_status_t status;
    status = mkl_sparse_d_create_csr(A, SPARSE_INDEX_BASE_ZERO, M.cols(), M.rows(), M.outerIndexPtr(), M.outerIndexPtr() + 1, M.innerIndexPtr(), M.valuePtr());
    if (status != 0)
    {
        std::cout << "CREATING MKL SPARSE MATRIX FAILED, CODE: " << status << "\n";
    }
}

IGL_INLINE void create_mkl_csr_matrix(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, sparse_matrix_t *A)
{
    sparse_status_t status;
    status = mkl_sparse_d_create_csr(A, SPARSE_INDEX_BASE_ZERO, M.rows(), M.cols(), M.outerIndexPtr(), M.outerIndexPtr() + 1, M.innerIndexPtr(), M.valuePtr());
    if (status != 0)
    {
        std::cout << "CREATING MKL SPARSE MATRIX FAILED, CODE: " << status << "\n";
    }
}
IGL_INLINE void mkl_export_csr(const sparse_matrix_t A,
                               sparse_index_base_t *indexing, MKL_INT *num_rows, MKL_INT *num_cols,
                               MKL_INT **row_start, MKL_INT **row_end, MKL_INT **col_indx,
                               double **values)
{
    sparse_status_t status;
    status = mkl_sparse_d_export_csr(A, indexing, num_rows, num_cols,
                                     row_start, row_end,
                                     col_indx, values);
    if (status != 0)
    {
        std::cout << "MKL EXPORT CSR FAILED, CODE: " << status << "\n";
    }
}

IGL_INLINE void export_csr_from_mkl(sparse_matrix_t &mkl_matrix, std::vector<int> &outerindex, std::vector<int> &innerindex, std::vector<double> &valueptr)
{
    sparse_status_t status;
    sparse_index_base_t indexing;
    MKL_INT num_rows, num_cols;
    MKL_INT *row_start, *row_end, *col_indx;
    double *values;
    mkl_export_csr(mkl_matrix, &indexing, &num_rows, &num_cols,
                            &row_start, &row_end,
                            &col_indx, &values);
    int count = 0;
    outerindex.resize(0);
    innerindex.resize(0);
    valueptr.resize(0);
    outerindex.reserve(num_rows + 1);
    innerindex.reserve(num_rows);
    valueptr.reserve(num_rows);
    for (int i = 0; i < num_rows; ++i)
    {
        outerindex.push_back(row_start[i]);
        std::vector<std::pair<int, double>> sorted_index;
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
}

IGL_INLINE void export_val_from_mkl(sparse_matrix_t &mkl_matrix, std::vector<double> &valueptr)
{
    sparse_status_t status;
    sparse_index_base_t indexing;
    MKL_INT num_rows, num_cols;
    MKL_INT *row_start, *row_end, *col_indx;
    double *values;
    mkl_export_csr(mkl_matrix, &indexing, &num_rows, &num_cols,
                   &row_start, &row_end,
                   &col_indx, &values);
    int count = 0;
    int index = 0;
    for (int i = 0; i < num_rows; ++i)
    {
        std::vector<std::pair<int, double>> sorted_index;
        for (int j = row_start[i]; j < row_end[i]; ++j)
        {
            sorted_index.push_back({col_indx[j], values[count]});
            ++count;
        }
        sort(sorted_index.begin(), sorted_index.end());
        for (auto p : sorted_index)
        {
            valueptr[index] = std::get<1>(p);
            index++;
        }
    }
}

IGL_INLINE void mkl_sypr_pre(sparse_matrix_t A, sparse_matrix_t B, struct matrix_descr &symmetric_type, sparse_matrix_t *result)
{
    // struct matrix_descr symmetric_type;
    symmetric_type.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    symmetric_type.mode = SPARSE_FILL_MODE_UPPER;
    symmetric_type.diag = SPARSE_DIAG_NON_UNIT;
    sparse_status_t status;
    status = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, A, B, symmetric_type, result, SPARSE_STAGE_NNZ_COUNT);
    if (status != 0)
    {
        std::cout << "SPARSE STAGE NNZ COUNT FAILED, CODE: " << status << "\n";
    }
    status = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, A, B, symmetric_type, result, SPARSE_STAGE_FINALIZE_MULT_NO_VAL);
    if (status != 0)
    {
        std::cout << "SPARSE STAGE FINALIZE MULT NO VAL FAILED, CODE: " << status << "\n";
    }
    mkl_set_num_threads_local(0);
}

IGL_INLINE void mkl_sypr_final(sparse_matrix_t A, sparse_matrix_t B, struct matrix_descr &symmetric_type, sparse_matrix_t *result)
{
    sparse_status_t status;
    mkl_set_num_threads_local(0);
    status = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, A, B, symmetric_type, result, SPARSE_STAGE_FINALIZE_MULT);
    if (status != 0)
    {
        std::cout << "SPARSE STAGE FINALIZE MULT FAILED, CODE: " << status << "\n";
    }
}

IGL_INLINE void mkl_sparse_add(sparse_matrix_t A, sparse_matrix_t B, double c, sparse_matrix_t *result)
{
    sparse_status_t status;
    mkl_set_num_threads_local(0);
    status = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, A, c, B, result);
    if (status != 0)
    {
        std::cout << "SPARSE ADD FAILED, CODE: " << status << "\n";
    }
}

} // namespace igl