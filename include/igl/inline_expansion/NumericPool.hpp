#include "NumericType.hpp"
#include <vector>
#include <map>
#pragma once

namespace ie {
class NumericPool {
public:
    NumericPool();
    static std::vector<NumericType> tree_node_pool;
    static std::map<double, size_t> const_node_position;
    void clear_pool();
};
}