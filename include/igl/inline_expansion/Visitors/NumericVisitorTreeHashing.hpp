#pragma once
#include "NumericVisitor.hpp"
#include "../NumericPool.hpp"
#include "../NumericType.hpp"
#include <map>
#include <mutex>
#include <vector>

namespace ie
{
class TreeToFileVisitor;
class TreeToIndexVisitor;
class NumericVisitorTreeHashing : public NumericVisitor
{
private:
    // a bool serves as a lock
    // needs to check with it to see if operation can be written at this moment
    static std::mutex can_write_operation_map;

    // another bool serves as a lock
    // needs to check with it to see if a constant can be written at this moment
    static std::mutex can_write_const_map;

    // an atomic operation to write to operation map
    size_t write_to_operation_map(const size_t &left_operation_id, const size_t &right_operation_id, const char &current_operation);

    // an atomic operation to write to constant map
    size_t write_to_constant_map(const double &const_value);

public:
    // default constructor
    NumericVisitorTreeHashing();

    // will record the operation type here, the key is in the format of:
    // <left operation key, right operation key>, current operation
    // current operation for leaf is 'i' and for constant is 'c'
    // the left oeration key for leaf will be the matrix number
    // for constant the left key is the seen constant position
    // right will be set to 0
    std::map<std::pair<std::pair<size_t, size_t>, char>, size_t> operation_to_id_map = {};
    std::vector<std::pair<std::pair<size_t, size_t>, char>> id_to_operation_map = {};

    // record for each tree in result array, what operation this element is
    std::vector<size_t> data_array_operation_ids;
    // record for each tree in result array, what data_ids this tree uses, note that we don't need
    // to record the matrix id because it will be stored in operation
    std::vector<std::vector<size_t>> data_array_used_data_ids;

    // storing the constants in it. the map is for looking up the position in seen_consts
    // because everything will be given a size_t as id
    std::map<double, size_t> seen_consts_position = {};
    std::vector<double> seen_consts = {};

    // visit function
    void visit(NumericType &n, size_t data_position, bool top_level, bool store_position) override;

    // accept function to generate file
    void accept(TreeToFileVisitor &visitor);

    // accept function to fill up indices
    void accept(TreeToIndexVisitor &visitor);

    // initialize with the length of the valuePtr of the result matrix
    NumericVisitorTreeHashing(unsigned int data_size);
};
} // namespace ie
