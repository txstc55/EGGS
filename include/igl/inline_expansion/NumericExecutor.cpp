#include "NumericExecutor.hpp"
#include <tbb/parallel_for.h>

using namespace std;

namespace ie
{
NumericExecutor::NumericExecutor() {}

NumericExecutor::NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::RowMajor> &R, size_t gap)
{
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(R.nonZeros());
    // tbb::parallel_for(size_t(0), size_t(R.nonZeros()), [&](size_t i) {
    //     R.valuePtr()[i].accept(this->th, i, true);
    // });
    for (int i = 0; i < R.nonZeros(); i++)
    {
        cout<<i<<" out of "<<R.nonZeros()<<"\n";
        R.valuePtr()[i].accept(this->th, i, true);
    }
    // this->choice = choice;
    // NumericType::clear_pool();

    // for grouped function
    cout << "Group by tree type method chosen\n";
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    this->tg.compile_file();
    this->tg.link_functions();
    cout << "Linking function completed\n";
}
NumericExecutor::NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor> &R, size_t gap)
{
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(R.nonZeros());
    // tbb::parallel_for(size_t(0), size_t(R.nonZeros()), [&](size_t i) {
    //     R.valuePtr()[i].accept(this->th, i, true);
    // });
    for (int i = 0; i < R.nonZeros(); i++)
    {
        R.valuePtr()[i].accept(this->th, i, true);
    }
    // this->choice = choice;
    // NumericType::clear_pool();

    // for grouped function
    cout << "Group by tree type method chosen\n";
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    this->tg.compile_file();
    this->tg.link_functions();
    cout << "Linking function completed\n";
}

NumericExecutor::NumericExecutor(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &R, size_t gap)
{
    int matrix_size = R.rows() * R.cols();
    this->gap = gap;
    this->th = NumericVisitorTreeHashing(matrix_size);
    // tbb::parallel_for(size_t(0), size_t(matrix_size), [&](size_t i) {
    //     R(i % R.rows(), i / R.rows()).accept(this->th, i, true);
    // });
    for (int i = 0; i < matrix_size; i++)
    {
        R(i % R.rows(), i / R.rows()).accept(this->th, i, true);
    }
    cout << "Group by tree type method chosen\n";
    this->tg = TreeToFileVisitorGroupByFunction(gap);
    this->th.accept(this->tg);
    this->tg.compile_file();
    this->tg.link_functions();
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