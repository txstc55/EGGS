#include "NumericExecutor.hpp"
#include <tbb/parallel_for.h>
#include <iostream>
#include "../Timer.h"
using namespace std;

namespace ie
{
NumericExecutor::NumericExecutor() {}

NumericExecutor::NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::RowMajor> &R, size_t gap)
{
    igl::Timer t;
    t.start();
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(R.nonZeros());
    for (int i = 0; i < R.nonZeros(); i++)
    {
        R.valuePtr()[i].accept(this->th, i, true);
    }
    // for grouped function
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    std::cout<<"Code generation took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    t.start();
    this->tg.compile_file();
    this->tg.link_functions();
    std::cout<<"Compiling code took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    cout << "Linking function completed\n";
}
NumericExecutor::NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> &R, size_t gap)
{
    igl::Timer t;
    t.start();
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(R.nonZeros());
    for (int i = 0; i < R.nonZeros(); i++)
    {
        R.valuePtr()[i].accept(this->th, i, true);
    }

    // for grouped function
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    std::cout<<"Code generation took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    t.start();
    this->tg.compile_file();
    this->tg.link_functions();
    std::cout<<"Compiling code took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    cout << "Linking function completed\n";
}

NumericExecutor::NumericExecutor(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &R, size_t gap)
{
    igl::Timer t;
    t.start();
    int matrix_size = R.rows() * R.cols();
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(matrix_size);
    for (int i = 0; i < matrix_size; i++)
    {
        R(i % R.rows(), i / R.rows()).accept(this->th, i, true);
    }
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    std::cout<<"Code generation took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    t.start();
    this->tg.compile_file();
    this->tg.link_functions();
    std::cout<<"Compiling code took: "<<t.getElapsedTimeInMilliSec()<<" ms\n";
    cout << "Linking function completed\n";
}

void NumericExecutor::ExecuteSingle(const vector<vector<double>> &data, vector<double> &result)
{
    this->tg.Execute_single(this->tg.reordered_data_id, this->tg.reordered_result_position, data, result);
}

void NumericExecutor::ExecuteMulti(const vector<vector<double>> &data, vector<double> &result)
{
    this->tg.Execute_multi(this->tg.reordered_data_id, this->tg.reordered_result_position, data, result);
}
} // namespace ie