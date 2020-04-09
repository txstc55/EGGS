#include "NumericType.hpp"
#include <vector>
#pragma once

namespace ie {
class NumericPool {
public:
    NumericPool();
    static std::vector<NumericType> tree_node_pool;
    void clear_pool();
};
}