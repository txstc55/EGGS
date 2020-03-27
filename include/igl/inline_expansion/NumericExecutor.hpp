#pragma once
#include "GroupByFunction/TreeToFileVisitorGroupByFunction.hpp"
#include "NumericType.hpp"
#include "Visitors/NumericVisitor.hpp"
#include <Eigen/Sparse>

namespace ie {
class NumericExecutor {
private:
    NumericVisitorTreeHashing th;

    // group by function, each time two element is executed
    TreeToFileVisitorGroupByFunction tg;

    size_t gap;
    size_t choice = 0;

public:
    NumericExecutor();
    NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::RowMajor>& R, size_t gap);
    NumericExecutor(Eigen::SparseMatrix<ie::NumericType, Eigen::ColMajor>& R, size_t gap);
    NumericExecutor(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic>& R, size_t gap);
    void ExecuteSingle(const std::vector<std::vector<double>>& data, std::vector<double>& result);
    void ExecuteMulti(const std::vector<std::vector<double>>& data, std::vector<double>& result);
};
}