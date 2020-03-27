#pragma once
#include "../Visitors/TreeToFileVisitor.hpp"
#include <set>

namespace ie
{
class TreeToFileVisitorGroupByFunction: public TreeToFileVisitor {
private:
    size_t gap = 0;


public:
    // init
    TreeToFileVisitorGroupByFunction();

    TreeToFileVisitorGroupByFunction(size_t gap);

    // the visit function
    void visit(NumericVisitorTreeHashing& trees) override;


    // the reordered data id usage
    std::vector<size_t> reordered_data_id;
    std::vector<size_t> reordered_result_position;

    void link_functions() override;

    void (*Execute_single)(const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<std::vector<double>>&, std::vector<double>&);
    void (*Execute_multi)(const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<std::vector<double>>&, std::vector<double>&);

    std::string function_hashed_name;

};




}