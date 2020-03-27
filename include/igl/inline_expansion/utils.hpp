#ifndef INLINE_EXPANSION_UTILS_HPP_
#define INLINE_EXPANSION_UTILS_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <unordered_set>

#include "NumericPool.hpp"
#include "NumericType.hpp"

#include <cstdlib>
#include <fstream>
#include <string>
#include <tbb/parallel_for.h>
#include <vector>

#include <iostream>
#include <tbb/task_scheduler_init.h>
#include <x86intrin.h>
#include <xmmintrin.h> // SSE

namespace ie {

/*!
     * Generate a dense matrix of size m x n
     */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
generate_dense_matrix(int m, int n)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = rand() / ((T)RAND_MAX);
        }
    }
    return A;
}

/*!
     * Generate a sparse matrix of size m x n, with density p
     */
template <typename T, int _Options = Eigen::RowMajor>
Eigen::SparseMatrix<T>
generate_sparse_matrix(int m, int n, int entry_per_row)
{
    Eigen::SparseMatrix<T, _Options> A(m, n);
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(m * entry_per_row);
    for (int i = 0; i < m; ++i) {
        std::set<int> col_nums;
        while (col_nums.size() < entry_per_row) {
            col_nums.insert(int(rand() % n));
        }

        for (auto c : col_nums) {
            trip.push_back(Eigen::Triplet<double>(i, c, rand() / ((T)RAND_MAX)));
        }

    }
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    return A;
}

/*!
     * Generate a sparse matrix of size m x n, with average density p
     */
template <typename T, int _Options = Eigen::RowMajor>
Eigen::SparseMatrix<T>
generate_sparse_matrix_average(int m, int n, int entry_per_row)
{
    Eigen::SparseMatrix<T, _Options> A(m, n);
    std::set<std::pair<int, int>> positions;
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(m * entry_per_row);
    while (positions.size() < entry_per_row * m) {
        positions.insert({ int(rand() % m), int(rand() % n) });
    }
    for (auto p : positions) {
        trip.push_back(Eigen::Triplet<double>(p.first, p.second, rand() / ((T)RAND_MAX)));
    }
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    return A;
}

// generate a sparse vector
template <typename T, int _Options = Eigen::RowMajor>
Eigen::SparseMatrix<T>
generate_sparse_vector(int m, double density)
{
    Eigen::SparseMatrix<T, _Options> A(m, 1);
    int number_of_entries = 0;
    if (density > 0.5) {
        number_of_entries = m * (1 - density);
    } else {
        number_of_entries = m * density;
    }
    std::set<int> positions;
    while (positions.size() < number_of_entries) {
        positions.insert(int(rand() % m));
    }
    if (density <= 0.5) {
        for (auto p : positions) {
            A.insert(p, 0) = rand() / ((T)RAND_MAX);
        }
    } else {
        for (int i = 0; i < m; i++) {
            if (positions.find(i) == positions.end()) {
                A.insert(i, 0) = rand() / ((T)RAND_MAX);
            }
        }
    }
    A.makeCompressed();
    return A;
}

// generate a symmetric matrix
template <typename T, int _Options = Eigen::RowMajor>
Eigen::SparseMatrix<T>
generate_sparse_symmetric_matrix(int m, int entry_per_row)
{
    Eigen::SparseMatrix<T, _Options> A(m, m);
    std::vector<std::vector<std::pair<int, double>>> row_index_value;
    std::vector<std::unordered_set<int>> row_unique_col;
    row_index_value.resize(m);
    row_unique_col.resize(m);
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(m * entry_per_row);
    for (int i = 0; i < m; ++i) {
        while (row_unique_col[i].size() < (entry_per_row / 2)) {
            int random_col = int(rand() % m);
            if (random_col != i && row_unique_col[random_col].find(i) == row_unique_col[random_col].end()) {
                row_unique_col[i].insert(random_col);
            }
        }
    }
    for (int i = 0; i < m; i++) {
        for (const int ele : row_unique_col[i]) {
            T random_value = rand() / ((T)RAND_MAX);
            row_index_value[i].push_back({ ele, random_value });
            row_index_value[ele].push_back({ i, random_value });
        }
    }
    int count = 0;
    for (int i = 0; i < m; i++) {
        for (auto it : row_index_value[i]) {
            trip.push_back(Eigen::Triplet<double>(i, it.first, it.second));
        }
        trip.push_back(Eigen::Triplet<double>(i, i, rand() / ((T)RAND_MAX)));
        count += row_index_value[i].size() + 1;
    }
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    return A;
}

/*!
     * \brief Convert a user-given sample sparse matrix to a NumericType matrix
     *        which has the same sparse structure and storage order (RowMajor or
     *        ColumnMajor) as original matrix
     * \param A dense matrix to be converted
     * \param matrix_id The index of M in the original function (0-indexed)
     */
template <typename T, int _Options = Eigen::RowMajor>
Eigen::SparseMatrix<ie::NumericType, _Options>
to_sparse_numeric(
    const Eigen::SparseMatrix<T, _Options>& A, int matrix_id)
{
    Eigen::SparseMatrix<ie::NumericType, _Options> M(A.rows(), A.cols());
    unsigned int data_id = 0;
    std::vector<Eigen::Triplet<ie::NumericType>> trip;
    trip.reserve(A.nonZeros());
    (*ie::NumericType::pool).tree_node_pool.reserve(std::max((*ie::NumericType::pool).tree_node_pool.size() + A.nonZeros(), (*ie::NumericType::pool).tree_node_pool.size() * 2 * int(A.nonZeros() / A.rows())));
    for (unsigned int i = 0; i < A.outerSize(); ++i) {
        for (typename Eigen::SparseMatrix<T, _Options>::InnerIterator it(A, i); it; ++it) {
            trip.push_back(Eigen::Triplet<ie::NumericType>(it.row(), it.col(), ie::NumericType(matrix_id, data_id)));
            ++data_id;
        }
    }
    M.setFromTriplets(trip.begin(), trip.end());
    M.makeCompressed();
    return M;
}

/*!
     * \brief Convert a user-given sample dense matrix to a NumericType matrix
     *        which has the same size and storage order (RowMajor or ColumnMajor)
     *        as original matrix
     * \param A dense matrix to be converted
     * \param matrix_id The index of M in the original function (0-indexed)
     */
template <typename T>
Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic>
to_dense_numeric(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
    int matrix_id)
{
    Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> M(A.rows(), A.cols());
    (*ie::NumericType::pool).tree_node_pool.reserve((*ie::NumericType::pool).tree_node_pool.size() + A.rows() * A.cols());
    int count = 0;
    for (int i = 0; i < A.cols(); ++i) {
        for (int j = 0; j < A.rows(); ++j) {
            // set matrix id and fake data id
            M(j, i) = ie::NumericType(matrix_id, count);
            count++;
        }
    }

    return M;
}

}

#endif
