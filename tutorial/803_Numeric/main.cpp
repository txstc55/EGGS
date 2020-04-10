////////////////////////////////////////////////////////////////////////////////
// Keep this file empty, and implement unit tests in separate compilation units!
////////////////////////////////////////////////////////////////////////////////

#include <array>

#include <igl/inline_expansion/utils.hpp>
#include "test_utils.hpp"
#include <string>
#include "mkl.h"
#include "mkl_spblas.h"
#include <fstream>
#include <iomanip>

#include <igl/inline_expansion/NumericExecutor.hpp>
#include "tbb/blocked_range.h"
#include <igl/Timer.h>
#include <fstream>

#include <algorithm>
#include <math.h>

using namespace std;
using namespace Eigen;
using namespace ie;

// ========================================================================================
// ========================================================================================
// ========================================================================================
template <class ValueType>
ValueType read_option(const char *option, int argc, char **argv, const char *default_value = nullptr);

template <>
std::string read_option<std::string>(const char *option, int argc, char **argv, const char *default_value)
{
    for (int i = 0; i < argc - 1; i++)
    {
        if (!strcmp(argv[i], option))
        {
            return std::string(argv[i + 1]);
        }
    }
    if (default_value)
        return std::string(default_value);
    std::cerr << "Option " << option << " was not provided. Exiting...\n";
    exit(1);
}
template <>
int read_option<int>(const char *option, int argc, char **argv, const char *default_value)
{
    return strtol(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL, 10);
}
template <>
long read_option<long>(const char *option, int argc, char **argv, const char *default_value)
{
    return strtol(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL, 10);
}
template <>
float read_option<float>(const char *option, int argc, char **argv, const char *default_value)
{
    return strtod(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL);
}
template <>
double read_option<double>(const char *option, int argc, char **argv, const char *default_value)
{
    return strtof(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL);
}
// ========================================================================================
// ========================================================================================
// ========================================================================================

int main(int argc, char *argv[])
{
    int matrix_size = read_option<int>("-r", argc, argv, "100000");
    int entry_per_row = read_option<int>("-e", argc, argv, "5");
    int demo_type = read_option<int>("-d", argc, argv, "0"); // 0 for A*B, 1 for operations involves constants, 2 for SYPR, 3 for SYRK
    sparse_matrix_t R1_mkl_t, R2_mkl_t;                      // store result for multi thread and single thread
    vector<double> numeric_result1, numeric_result2;         // to store the numeric test results for multi thread and single thread
    NumericExecutor ex;                                      // the executor
    srand(1);                                                // fix the seed
    ofstream result_file;                                    // the output file for time result
    switch (demo_type)
    {
    case 0:
    {
        result_file.open("pipe_result.txt", std::ios_base::app);
        SparseMatrix<double, RowMajor> m1 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        SparseMatrix<double, RowMajor> m2 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        cout << "Matrix generated for this operation\n";
        igl::Timer numeric_prep;
        numeric_prep.start();
        SparseMatrix<NumericType, RowMajor> m1_numeric = to_sparse_numeric<double, RowMajor>(m1, 0);
        SparseMatrix<NumericType, RowMajor> m2_numeric = to_sparse_numeric<double, RowMajor>(m2, 1);
        SparseMatrix<NumericType, RowMajor> result_numeric = m1_numeric * m2_numeric;
        ex = NumericExecutor(result_numeric, 0);
        numeric_prep.stop();
        cout << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        result_file << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        vector<vector<double>> DATAS = {vector<double>(m1.valuePtr(), m1.valuePtr() + m1.nonZeros()), vector<double>(m2.valuePtr(), m2.valuePtr() + m2.nonZeros())};
        vector<SparseMatrix<double, RowMajor>> MATRIX_VECTOR = {m1, m2};
        PROFILE_MKL_MULTI_PIPE(MATRIX_VECTOR, R1_mkl_t, result_file);
        PROFILE_MKL_SINGLE_PIPE(MATRIX_VECTOR, R2_mkl_t, result_file);
        numeric_result1.resize(result_numeric.nonZeros());
        numeric_result2.resize(result_numeric.nonZeros());
        PROFILE_EXECUTOR_MULTI(ex, DATAS, numeric_result1, result_file);
        PROFILE_EXECUTOR_SINGLE(ex, DATAS, numeric_result2, result_file);
        // extract the result
        SparseMatrix<double, RowMajor> Eigen_result;
        auto elapsed_eigen_single = benchmarkTimer([&]() {
            for (int i = 0; i < 100; i++)
            {
                Eigen_result = m1 * m2;
            }
        });
        cout << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        result_file << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        auto R1_info = extract_value<double>(R1_mkl_t);
        auto R_mkl = ConstructSparseMatrix(get<0>(R1_info), get<1>(R1_info), get<2>(R1_info), (get<5>(R1_info)).data(), (get<3>(R1_info)).data(), (get<4>(R1_info)).data());
        auto R_numeric = ConstructSparseMatrix(result_numeric.rows(), result_numeric.cols(), result_numeric.nonZeros(), numeric_result2.data(), result_numeric.outerIndexPtr(), result_numeric.innerIndexPtr());
        cout << "MKL ERROR: " << (Eigen_result - R_mkl).norm() << "\n";
        cout << "NUMERIC ERROR: " << (Eigen_result - R_numeric).norm() << "\n";
        break;
    }
    case 1:
    {
        result_file.open("sypr_result.txt", std::ios_base::app);
        SparseMatrix<double, RowMajor> m1 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        SparseMatrix<double, RowMajor> m2 = generate_sparse_symmetric_matrix<double, RowMajor>(matrix_size, 1);
        cout << "Matrix generated for this operation\n";
        igl::Timer numeric_prep;
        numeric_prep.start();
        SparseMatrix<NumericType, RowMajor> m1_numeric = to_sparse_numeric<double, RowMajor>(m1, 0);
        SparseMatrix<NumericType, RowMajor> m2_numeric = to_sparse_numeric<double, RowMajor>(m2, 1);
        SparseMatrix<NumericType, RowMajor> result_numeric = (SparseMatrix<NumericType, RowMajor>(m1_numeric.transpose()) * m2_numeric * m1_numeric).triangularView<Upper>();
        ex = NumericExecutor(result_numeric, 0);
        numeric_prep.stop();
        cout << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        result_file << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        vector<vector<double>> DATAS = {vector<double>(m1.valuePtr(), m1.valuePtr() + m1.nonZeros()), vector<double>(m2.valuePtr(), m2.valuePtr() + m2.nonZeros())};
        vector<SparseMatrix<double, RowMajor>> MATRIX_VECTOR = {m1, m2};
        PROFILE_MKL_MULTI_SYPR(MATRIX_VECTOR, R1_mkl_t, result_file);
        PROFILE_MKL_SINGLE_SYPR(MATRIX_VECTOR, R2_mkl_t, result_file);
        numeric_result1.resize(result_numeric.nonZeros());
        numeric_result2.resize(result_numeric.nonZeros());
        PROFILE_EXECUTOR_MULTI(ex, DATAS, numeric_result1, result_file);
        PROFILE_EXECUTOR_SINGLE(ex, DATAS, numeric_result2, result_file);
        SparseMatrix<double, RowMajor> Eigen_result;
        auto elapsed_eigen_single = benchmarkTimer([&]() {
            for (int i = 0; i < 100; i++)
            {
                Eigen_result = SparseMatrix<double, RowMajor>(m1.transpose()) * m2 * m1;
            }
        });
        cout << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        result_file << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        auto R1_info = extract_value<double>(R1_mkl_t);
        auto R_mkl = ConstructSparseMatrix(get<0>(R1_info), get<1>(R1_info), get<2>(R1_info), (get<5>(R1_info)).data(), (get<3>(R1_info)).data(), (get<4>(R1_info)).data());
        auto R_numeric = ConstructSparseMatrix(result_numeric.rows(), result_numeric.cols(), result_numeric.nonZeros(), numeric_result2.data(), result_numeric.outerIndexPtr(), result_numeric.innerIndexPtr());
        cout << "MKL ERROR: " << (Eigen_result.triangularView<Upper>() - R_mkl).norm() << "\n";
        cout << "NUMERIC ERROR: " << (Eigen_result.triangularView<Upper>() - R_numeric).norm() << "\n";
        break;
    }
    case 2:
    {
        result_file.open("syrk_result.txt", std::ios_base::app);
        SparseMatrix<double, RowMajor> m1 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        cout << "Matrix generated for this operation\n";
        igl::Timer numeric_prep;
        numeric_prep.start();
        SparseMatrix<NumericType, RowMajor> m1_numeric = to_sparse_numeric<double, RowMajor>(m1, 0);
        SparseMatrix<NumericType, RowMajor> result_numeric = (SparseMatrix<NumericType, RowMajor>(m1_numeric.transpose()) * m1_numeric).triangularView<Upper>();
        ex = NumericExecutor(result_numeric, 0);
        numeric_prep.stop();
        cout << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        result_file << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        vector<vector<double>> DATAS = {vector<double>(m1.valuePtr(), m1.valuePtr() + m1.nonZeros())};
        vector<SparseMatrix<double, RowMajor>> MATRIX_VECTOR = {m1};
        PROFILE_MKL_MULTI_SYRK(MATRIX_VECTOR, R1_mkl_t, result_file);
        PROFILE_MKL_SINGLE_SYRK(MATRIX_VECTOR, R2_mkl_t, result_file);
        numeric_result1.resize(result_numeric.nonZeros());
        numeric_result2.resize(result_numeric.nonZeros());
        PROFILE_EXECUTOR_MULTI(ex, DATAS, numeric_result1, result_file);
        PROFILE_EXECUTOR_SINGLE(ex, DATAS, numeric_result2, result_file);
        SparseMatrix<double, RowMajor> Eigen_result;
        auto elapsed_eigen_single = benchmarkTimer([&]() {
            for (int i = 0; i < 100; i++)
            {
                Eigen_result = SparseMatrix<double, RowMajor>(m1.transpose()) * m1;
            }
        });
        cout << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        result_file << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        auto R1_info = extract_value<double>(R1_mkl_t);
        auto R_mkl = ConstructSparseMatrix(get<0>(R1_info), get<1>(R1_info), get<2>(R1_info), (get<5>(R1_info)).data(), (get<3>(R1_info)).data(), (get<4>(R1_info)).data());
        auto R_numeric = ConstructSparseMatrix(result_numeric.rows(), result_numeric.cols(), result_numeric.nonZeros(), numeric_result2.data(), result_numeric.outerIndexPtr(), result_numeric.innerIndexPtr());
        cout << "MKL ERROR: " << (Eigen_result.triangularView<Upper>() - R_mkl).norm() << "\n";
        cout << "NUMERIC ERROR: " << (Eigen_result.triangularView<Upper>() - R_numeric).norm() << "\n";
        break;
    }
    case 3:
    {
        result_file.open("spmv_result.txt", std::ios_base::app);
        SparseMatrix<double, RowMajor> m1 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        Matrix<double, Dynamic, Dynamic> DENSE_VECTOR = generate_dense_matrix<double>(matrix_size, 1);
        cout << "Matrix generated for this operation\n";
        igl::Timer numeric_prep;
        numeric_prep.start();
        SparseMatrix<NumericType, RowMajor> m1_numeric = to_sparse_numeric<double, RowMajor>(m1, 0);
        Matrix<NumericType, Dynamic, Dynamic> DENSE_VECTOR_numeric = to_dense_numeric(DENSE_VECTOR, 1);
        cout << "Converting finished\n";
        Matrix<NumericType, Dynamic, Dynamic> result_numeric = m1_numeric * DENSE_VECTOR_numeric;
        ex = NumericExecutor(result_numeric, 0);
        numeric_prep.stop();
        cout << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        result_file << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        vector<vector<double>> DATAS = {vector<double>(m1.valuePtr(), m1.valuePtr() + m1.nonZeros()), vector<double>(DENSE_VECTOR.data(), DENSE_VECTOR.data() + DENSE_VECTOR.rows())};
        Matrix<double, Dynamic, 1> mkl_result;
        mkl_result.resize(matrix_size, 1);
        PROFILE_MKL_MULTI_SPMV(m1, DENSE_VECTOR.data(), mkl_result.data(), result_file);
        numeric_result1.resize(matrix_size);
        PROFILE_EXECUTOR_MULTI(ex, DATAS, numeric_result1, result_file);
        MatrixXd numeric_result;
        SparseMatrix<double, RowMajor> Eigen_result;
        auto elapsed_eigen_single = benchmarkTimer([&]() {
            for (int i = 0; i < 100; i++)
            {
                Eigen_result = SparseMatrix<double, RowMajor>(m1.transpose()) * m1;
            }
        });
        cout << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        result_file << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        numeric_result.resize(matrix_size, 1);
        numeric_result = Map<MatrixXd>(numeric_result1.data(), matrix_size, 1);
        Matrix<double, Dynamic, Dynamic> Eigen_result;
        cout << "MKL ERROR: " << (Eigen_result - mkl_result).norm() << "\n";
        cout << "NUMERIC ERROR: " << (Eigen_result - numeric_result).norm() << "\n";
        break;
    }
    case 4:
    {
        result_file.open("const_result.txt", std::ios_base::app);
        SparseMatrix<double, RowMajor> m1 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        SparseMatrix<double, RowMajor> m2 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        SparseMatrix<double, RowMajor> m3 = generate_sparse_matrix_average<double, RowMajor>(matrix_size, matrix_size, entry_per_row);
        cout << "Matrix generated for this operation\n";
        igl::Timer numeric_prep;
        numeric_prep.start();
        SparseMatrix<NumericType, RowMajor> m1_numeric = to_sparse_numeric<double, RowMajor>(m1, 0);
        SparseMatrix<NumericType, RowMajor> m2_numeric = to_sparse_numeric<double, RowMajor>(m2, 1);
        SparseMatrix<NumericType, RowMajor> m3_numeric = to_sparse_numeric<double, RowMajor>(m3, 2);
        SparseMatrix<NumericType, RowMajor> result_numeric = SparseMatrix<NumericType, RowMajor>((3.78 * m1_numeric + m2_numeric).transpose()) * (6.942 * SparseMatrix<NumericType, RowMajor>(m2_numeric.transpose()) + m3_numeric);
        ex = NumericExecutor(result_numeric, 0);
        numeric_prep.stop();
        cout << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        result_file << "Numeric pre-computation: " << numeric_prep.getElapsedTimeInMicroSec() << " us\n";
        vector<vector<double>> DATAS = {vector<double>(m1.valuePtr(), m1.valuePtr() + m1.nonZeros()), vector<double>(m2.valuePtr(), m2.valuePtr() + m2.nonZeros()), vector<double>(m3.valuePtr(), m3.valuePtr() + m3.nonZeros())};
        vector<SparseMatrix<double, RowMajor>> MATRIX_VECTOR = {m1, m2, m3};
        PROFILE_MKL_CONST_OP_MULT(MATRIX_VECTOR, R1_mkl_t, result_file);
        PROFILE_MKL_CONST_OP_SINGLE(MATRIX_VECTOR, R2_mkl_t, result_file);
        numeric_result1.resize(result_numeric.nonZeros());
        numeric_result2.resize(result_numeric.nonZeros());
        PROFILE_EXECUTOR_MULTI(ex, DATAS, numeric_result1, result_file);
        PROFILE_EXECUTOR_SINGLE(ex, DATAS, numeric_result2, result_file);
        // extract the result
        SparseMatrix<double, RowMajor> Eigen_result;
        auto elapsed_eigen_single = benchmarkTimer([&]() {
            for (int i = 0; i < 100; i++)
            {
                Eigen_result = SparseMatrix<double, RowMajor>((3.78 * m1 + m2).transpose()) * (6.942 * SparseMatrix<double, RowMajor>(m2.transpose()) + m3);
            }
        });
        cout << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        result_file << "EIGEN SINGLE THREAD: " << elapsed_eigen_single << " us\n";
        auto R1_info = extract_value<double>(R1_mkl_t);
        auto R_mkl = ConstructSparseMatrix(get<0>(R1_info), get<1>(R1_info), get<2>(R1_info), (get<5>(R1_info)).data(), (get<3>(R1_info)).data(), (get<4>(R1_info)).data());
        auto R_numeric = ConstructSparseMatrix(result_numeric.rows(), result_numeric.cols(), result_numeric.nonZeros(), numeric_result2.data(), result_numeric.outerIndexPtr(), result_numeric.innerIndexPtr());
        cout << "MKL ERROR: " << (Eigen_result - R_mkl).norm() << "\n";
        cout << "NUMERIC ERROR: " << (Eigen_result - R_numeric).norm() << "\n";
        break;
    }
    }
    result_file.close();
}
