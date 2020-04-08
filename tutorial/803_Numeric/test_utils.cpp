#include "test_utils.hpp"

using namespace Eigen;
using namespace std;


template<>
sparse_status_t mkl_create_coo<float>(sparse_matrix_t* A,
                                      sparse_index_base_t indexing, MKL_INT num_rows, MKL_INT num_cols,
                                      MKL_INT nnz, MKL_INT* rows, MKL_INT* cols, float* values) {
    return mkl_sparse_s_create_coo(A, indexing, num_rows, num_cols,
                                   nnz, rows, cols, values);
}

template<>
sparse_status_t mkl_create_coo<double>(sparse_matrix_t* A,
                                       sparse_index_base_t indexing, MKL_INT num_rows, MKL_INT num_cols,
                                       MKL_INT nnz, MKL_INT* rows, MKL_INT* cols, double* values) {
    return mkl_sparse_d_create_coo(A, indexing, num_rows, num_cols,
                                   nnz, rows, cols, values);
}

template<>
sparse_status_t mkl_export_csr<float>(const sparse_matrix_t A,
                                      sparse_index_base_t* indexing, MKL_INT* num_rows, MKL_INT* num_cols,
                                      MKL_INT** row_start, MKL_INT** row_end, MKL_INT** col_indx,
                                      float** values) {
    return mkl_sparse_s_export_csr(A, indexing, num_rows, num_cols,
                                   row_start, row_end,
                                   col_indx, values);
}

template<>
sparse_status_t mkl_export_csr<double>(const sparse_matrix_t A,
                                       sparse_index_base_t* indexing, MKL_INT* num_rows, MKL_INT* num_cols,
                                       MKL_INT** row_start, MKL_INT** row_end, MKL_INT** col_indx,
                                       double** values) {
    return mkl_sparse_d_export_csr(A, indexing, num_rows, num_cols,
                                   row_start, row_end,
                                   col_indx, values);
}

template<>
void create_mkl_csr_matrix(SparseMatrix<double, RowMajor>& M, sparse_matrix_t* A) {
    sparse_status_t status;
    status = mkl_sparse_d_create_csr(A, SPARSE_INDEX_BASE_ZERO, M.rows(), M.cols(), M.outerIndexPtr(), M.outerIndexPtr() + 1, M.innerIndexPtr(), M.valuePtr());
}

template<>
void create_mkl_csr_matrix(SparseMatrix<float, RowMajor>& M, sparse_matrix_t* A) {
    sparse_status_t status;
    status = mkl_sparse_s_create_csr(A, SPARSE_INDEX_BASE_ZERO, M.rows(), M.cols(), M.outerIndexPtr(), M.outerIndexPtr() + 1, M.innerIndexPtr(), M.valuePtr());

}






