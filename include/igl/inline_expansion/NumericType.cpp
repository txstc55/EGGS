#include "NumericType.hpp"
#include "NumericPool.hpp"
#include "Visitors/NumericVisitor.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;

namespace ie
{
// vector<tuple<bool, bool, unsigned int, unsigned int, char>> NumericType::memory_pool = {};
// map<double, int> NumericType::existed_const_positions = {};
// vector<double> NumericType::const_pool = {};
NumericPool *NumericType::pool;

map<NumericType::NodeType, char> NumericType::node_to_op = {{NumericType::Constant, 'c'},
                                                            {NumericType::Leaf, 'i'},
                                                            {NumericType::Add, '+'},
                                                            {NumericType::Subtract, '-'},
                                                            {NumericType::Divide, '/'},
                                                            {NumericType::Multiply, '*'},
                                                            {NumericType::Sqrt, 's'}};

char NumericType::char_for_operation()
{
    return NumericType::node_to_op[operation];
}

NumericType::NumericType(size_t m_id, size_t d_id)
{
    this->matrix_id = m_id;
    this->data_id = d_id;
    this->operation = Leaf;
    this->self_index = NumericType::pool->tree_node_pool.size();
    NumericType::pool->tree_node_pool.push_back(*this);
}

NumericType::NumericType(double v)
{
    this->operation = Constant;
    this->const_value = v;
    this->self_index = NumericType::pool->tree_node_pool.size();
    NumericType::pool->tree_node_pool.push_back(*this);
}

NumericType::NumericType(const NumericType &n)
{
    this->matrix_id = n.matrix_id;
    this->data_id = n.data_id;
    this->self_index = n.self_index;
    this->const_value = n.const_value;
    this->left_index = n.left_index;
    this->right_index = n.right_index;
    this->operation = n.operation;
}

NumericType::NumericType()
{
}

void NumericType::set_pool(NumericPool *p)
{
    NumericType::pool = p;
}

void NumericType::accept(NumericVisitor &nv, size_t data_position)
{
    nv.visit(*this, data_position);
}

// the operations will avoid using excessive amount of prentices
// left is this and right is v, always
NumericType NumericType::operator+(const NumericType &v) const
{
    if (this->operation == Constant && this->const_value == 0)
    {
        return v;
    }
    if (v.operation == Constant && v.const_value == 0)
    {
        return (*this);
    }
    NumericType parent = NumericType();
    parent.self_index = NumericType::pool->tree_node_pool.size();
    parent.left_index = this->self_index;
    parent.right_index = v.self_index;
    parent.operation = Add;
    NumericType::pool->tree_node_pool.push_back(parent);
    return parent;
}

NumericType NumericType::operator-(const NumericType &v) const
{
    if (v.operation == Constant && v.const_value == 0)
    {
        return (*this);
    }
    NumericType parent = NumericType();
    parent.self_index = NumericType::pool->tree_node_pool.size();
    parent.left_index = this->self_index;
    parent.right_index = v.self_index;
    parent.operation = Subtract;
    NumericType::pool->tree_node_pool.push_back(parent);
    return parent;
}

NumericType NumericType::operator*(const NumericType &v) const
{
    if (this->operation == Constant && this->const_value == 1)
    {
        return v;
    }
    if (v.operation == Constant && v.const_value == 1)
    {
        return (*this);
    }
    NumericType parent = NumericType();
    parent.self_index = NumericType::pool->tree_node_pool.size();
    parent.left_index = this->self_index;
    parent.right_index = v.self_index;
    parent.operation = Multiply;
    NumericType::pool->tree_node_pool.push_back(parent);
    return parent;
}

NumericType NumericType::operator/(const NumericType &v) const
{
    if (v.operation == Constant && v.const_value == 1)
    {
        return (*this);
    }
    NumericType parent = NumericType();
    parent.self_index = NumericType::pool->tree_node_pool.size();
    parent.left_index = this->self_index;
    parent.right_index = v.self_index;
    parent.operation = Divide;
    NumericType::pool->tree_node_pool.push_back(parent);
    return parent;
}

NumericType NumericType::operator+=(const NumericType &v)
{
    return (*this) = (*this) + v;
}

NumericType NumericType::sqrt() const
{
    NumericType parent = NumericType();
    parent.self_index = NumericType::pool->tree_node_pool.size();
    parent.left_index = this->self_index;
    parent.right_index = 0;
    parent.operation = Sqrt;
    NumericType::pool->tree_node_pool.push_back(parent);
    return parent;
}

NumericType NumericType::operator*(const double v)
{
    if (v == 1.0)
    {
        return (*this);
    }
    else
    {
        NumericType constant = NumericType(v);
        return ((*this) * constant);
    }
}

NumericType NumericType::operator/(const double v)
{
    if (v == 1.0)
    {
        return (*this);
    }
    else
    {
        NumericType constant = NumericType(v);
        return ((*this) / constant);
    }
}

NumericType operator*(const double v, const NumericType &n)
{
    if (v == 1.0)
    {
        return (n);
    }
    else
    {
        NumericType constant = NumericType(v);
        return (constant * n);
    }
}

bool NumericType::operator<=(const NumericType &v) const
{
    return false;
}

bool NumericType::operator<(const NumericType &v) const
{
    return false;
}

NumericType NumericType::change_matrix_id(size_t i)
{
    this->matrix_id = i;
    NumericType::pool->tree_node_pool[self_index].matrix_id = i;
    return *this;
}

NumericType NumericType::change_data_id(size_t i)
{
    this->data_id = i;
    NumericType::pool->tree_node_pool[self_index].data_id = i;
    return *this;
}

void NumericType::clear_pool()
{
    NumericType::pool->clear_pool();
}

} // namespace ie
