#pragma once
#include "mkl_types.h"
#include "mkl_spblas.h"

#include "igl_inline.h"

#include <Eigen/Sparse>

namespace igl
{

IGL_INLINE void create_mkl_csr_matrix(Eigen::SparseMatrix<double, Eigen::ColMajor> &M, sparse_matrix_t *A); // when the matrix is column major, you are just creating the transpose
IGL_INLINE void create_mkl_csr_matrix(Eigen::SparseMatrix<double, Eigen::RowMajor> &M, sparse_matrix_t *A);
IGL_INLINE void mkl_export_csr(const sparse_matrix_t A,
                               sparse_index_base_t *indexing, MKL_INT *num_rows, MKL_INT *num_cols,
                               MKL_INT **row_start, MKL_INT **row_end, MKL_INT **col_indx,
                               double **values);
IGL_INLINE void export_csr_from_mkl(sparse_matrix_t &mkl_matrix, std::vector<int> &outerindex, std::vector<int> &innerindex, std::vector<double> &valueptr);
IGL_INLINE void export_val_from_mkl(sparse_matrix_t &mkl_matrix, std::vector<double> &valueptr);
IGL_INLINE void mkl_sypr_pre(sparse_matrix_t A, sparse_matrix_t B, struct matrix_descr &symmetric_type, sparse_matrix_t *result);
IGL_INLINE void mkl_sypr_final(sparse_matrix_t A, sparse_matrix_t B, struct matrix_descr &symmetric_type, sparse_matrix_t *result);
IGL_INLINE void mkl_sparse_add(sparse_matrix_t A, sparse_matrix_t B, double c, sparse_matrix_t *result);

} // namespace igl