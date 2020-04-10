#include "NumericType.hpp"
#include "NumericPool.hpp"
#include "Visitors/NumericVisitor.hpp"
#include <cstdlib>
#include <iostream>
#include <queue>
#include <algorithm>
#include <set>

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
                                                            {NumericType::Sqrt, 's'},
                                                            {NumericType::Repeated, 'r'}};

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

void NumericType::accept(NumericVisitor &nv, size_t data_position, const bool top_level)
{
    map<size_t, size_t> chosen_repeated_node_map;
    (*this).MarkRepeatedNodes(chosen_repeated_node_map, top_level);       // we need to find the repeated node id first
    nv.visit(*this, data_position, true, true, chosen_repeated_node_map); // pass that to the visitor and do whatever they want to do with it
}

// the operations will avoid using excessive amount of prentices
// left is this and right is v, always
NumericType NumericType::operator+(const NumericType &v) const
{
    if (this->operation == Constant && this->const_value == 0)
        return v;
    if (v.operation == Constant && v.const_value == 0)
        return (*this);
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
        return (*this);
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
        return v;
    if (v.operation == Constant && v.const_value == 1)
        return (*this);
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
        return (*this);
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

NumericType NumericType::operator+(const double v)
{
    if (v == 0.0)
        return (*this);
    else
    {
        NumericType constant = NumericType(v);
        return ((*this) + constant);
    }
}

NumericType NumericType::operator*(const double v)
{
    if (v == 1.0)
        return (*this);
    else
    {
        NumericType constant = NumericType(v);
        return ((*this) * constant);
    }
}

NumericType NumericType::operator/(const double v)
{
    if (v == 1.0)
        return (*this);
    else
    {
        NumericType constant = NumericType(v);
        return ((*this) / constant);
    }
}

NumericType operator*(const double v, const NumericType &n)
{
    if (v == 1.0)
        return (n);
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

void NumericType::MarkRepeatedNodes(std::map<size_t, size_t> &chosen_repeated_node_map, const bool top_level)
{
    if (!top_level)
        return;
    queue<NumericType> candidates;
    vector<size_t> used_node_index;
    candidates.push(*this);
    // accumulate all the node index used
    // note that we do not care about repeated leaf index
    // and constant index
    // because those things can be repeated no matter what
    while (!candidates.empty())
    {
        NumericType c = candidates.front();
        candidates.pop();
        size_t left, right;
        switch (c.operation)
        {
        case Leaf:
            break;
        case Constant:
            break;
        case Sqrt:
        {
            left = c.left_index;
            candidates.push(NumericType::pool->tree_node_pool[left]);
            used_node_index.push_back(c.self_index);
            break;
        }
        case Add:
        case Subtract:
        case Divide:
        case Multiply:
        {
            left = c.left_index;
            right = c.right_index;
            candidates.push(NumericType::pool->tree_node_pool[left]);
            candidates.push(NumericType::pool->tree_node_pool[right]);
            used_node_index.push_back(c.self_index);
            break;
        }
        default:
            break;
        }
    }

    bool has_repeated = false;
    vector<size_t> repeated_indices;          // record what indices are repeated
    map<size_t, size_t> repeated_index_count; // record the actual repeated indices and how many times they appeared
    repeated_indices.reserve(used_node_index.size());
    for (unsigned int i = 0; i < used_node_index.size(); i++)
    {
        if (repeated_index_count.find(used_node_index[i]) == repeated_index_count.end())
            repeated_index_count.insert({used_node_index[i], 1});
        else
            repeated_index_count[used_node_index[i]]++;
    }
    for (auto p : repeated_index_count)
    {
        if (p.second > 1)
        {
            has_repeated = true;
            repeated_indices.push_back(p.first);
        }
    }
    // we have not detected any repeated entries
    if (!has_repeated)
        return;

    sort(repeated_indices.begin(), repeated_indices.end()); // sort the repeated indices

    std::set<size_t> chosen_repeated_node;
    set<size_t> chosen_repeated_node_children; // the choldren of the repeated nodes are repeated, but we do not want to label them
    // we have detected repeated entries, now do the work of finding the top level node that covers all of its children indices
    // printf("count %d \n", repeated_indices.size());
    for (int i = repeated_indices.size() - 1; i >= 0; i--)
    {
        // printf("at node %d\n", repeated_indices[i]);
        NumericType c = NumericType::pool->tree_node_pool[repeated_indices[i]];
        // we are not children of anyone
        if (chosen_repeated_node_children.find(repeated_indices[i]) == chosen_repeated_node_children.end())
        {
            size_t left, right;
            left = c.left_index;         // get the node used for left branch
            right = c.right_index;       // get the node used for right branch
            bool left_fulfilled = true;  // check if for n parents, we also have n left nodes
            bool right_fulfilled = true; // check if for n parents, we also have n right nodes
            if (left != right)
            {
                if (repeated_index_count.find(left) != repeated_index_count.end())
                {
                    if (repeated_index_count[left] > repeated_index_count[c.self_index])
                        left_fulfilled = false;
                }
                if (left_fulfilled && c.operation != Sqrt && repeated_index_count.find(right) != repeated_index_count.end())
                {
                    if (repeated_index_count[right] > repeated_index_count[c.self_index])
                        right_fulfilled = false;
                }
            }
            else
            {
                if (repeated_index_count.find(left) != repeated_index_count.end())
                {
                    if (repeated_index_count[left] > 2 * repeated_index_count[c.self_index])
                        left_fulfilled = false;
                }
            }
            // the occurrence of this node is the same
            if (left_fulfilled && right_fulfilled)
            {
                chosen_repeated_node.insert(c.self_index);          // this node can cover all the repeated nodes descended from itself
                chosen_repeated_node_children.insert(c.left_index); // we don't care that its children are repeated
                if (c.operation != Sqrt)
                    chosen_repeated_node_children.insert(c.right_index);
            }
        }
        else
        {
            chosen_repeated_node_children.insert(c.left_index); // insert the children of children so we don't end up putting grand children into repeated operation
            if (c.operation != Sqrt)
                chosen_repeated_node_children.insert(c.right_index);
        }
    }
    if (chosen_repeated_node.size() == 0)
        return;

    // now record from a bfs perspective, what order do we see those chosen nodes
    queue<NumericType> nodes;
    nodes.push(*this);
    while ((!nodes.empty()) && chosen_repeated_node.size() != chosen_repeated_node_map.size())
    {
        NumericType c = nodes.front();
        nodes.pop();
        int left, right;
        // check if it is in the chosen repeated nodes and if we have seen it
        switch (c.operation)
        {
        case Leaf:
            break;
        case Constant:
            break;
        case Sqrt:
        {
            if (chosen_repeated_node.find(c.self_index) != chosen_repeated_node.end())
            {
                if (chosen_repeated_node_map.find(c.self_index) == chosen_repeated_node_map.end())
                    chosen_repeated_node_map.insert({c.self_index, chosen_repeated_node_map.size()});
            }
            else
            {
                left = c.left_index;
                nodes.push(NumericType::pool->tree_node_pool[left]);
            }
            break;
        }
        case Add:
        case Subtract:
        case Divide:
        case Multiply:
        {
            if (chosen_repeated_node.find(c.self_index) != chosen_repeated_node.end())
            {
                if (chosen_repeated_node_map.find(c.self_index) == chosen_repeated_node_map.end())
                    chosen_repeated_node_map.insert({c.self_index, chosen_repeated_node_map.size()});
            }
            else
            {
                left = c.left_index;
                right = c.right_index;
                nodes.push(NumericType::pool->tree_node_pool[left]);
                nodes.push(NumericType::pool->tree_node_pool[right]);
            }
            break;
        }
        default:
            break;
        }
    }
    // printf("Number of chosen repeated nodes: %d\n", chosen_repeated_node_map.size());
}

} // namespace ie
