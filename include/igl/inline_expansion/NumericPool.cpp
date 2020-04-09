#include "NumericPool.hpp"
using namespace std;
namespace ie {
vector<NumericType> NumericPool::tree_node_pool = {};
map<double, size_t> NumericPool::const_node_position = {};

NumericPool::NumericPool()
{
    NumericType::set_pool(this);
}

void NumericPool::clear_pool() {
    NumericPool::tree_node_pool.clear();
    NumericPool::const_node_position.clear();
    NumericPool::tree_node_pool.shrink_to_fit();
}

}