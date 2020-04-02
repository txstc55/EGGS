#pragma once
#include <Eigen/Sparse>
#include <array>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <set>

namespace ie
{

class NumericPool;

class NumericVisitor;

class NumericType
{
private:
    // put everything into numeric type
    // get rid of the intermediate level
    // use inja to generate the file
public:
    // The type of the node
    enum NodeType
    {
        Leaf,
        Constant,
        Add,
        Subtract,
        Divide,
        Multiply,
        Sqrt,
        Repeated
    };

    // Return a char representation of the operation.
    // This should likely be removed in the future.
    char char_for_operation();

    // properties for a constant value
    double const_value;

    // record which matrix the data comes from
    size_t matrix_id;
    // record where the data comes from
    size_t data_id;

    // the type of operation
    NodeType operation;

    NumericType(size_t matrix_id, size_t data_id);
    NumericType(double v);
    NumericType(const NumericType &n);
    NumericType();

    // the pool for numeric type
    static NumericPool *pool;
    // set the pool pointer to something
    static void set_pool(NumericPool *p);
    // indices in the pool
    size_t self_index;
    size_t left_index;
    size_t right_index;

    // accept the visitor
    void accept(NumericVisitor &nv, size_t data_position, const bool top_level);

    // the operations will all return a NumericType
    // which has the structure of a tree
    NumericType operator+(const NumericType &v) const;
    NumericType operator-(const NumericType &v) const;
    NumericType operator*(const NumericType &v) const;
    NumericType operator/(const NumericType &v) const;
    NumericType operator+=(const NumericType &v);
    NumericType sqrt() const;

    NumericType operator*(const double v);
    NumericType operator/(const double v);

    friend NumericType operator*(const double v, const NumericType &n);

    bool operator<=(const NumericType &v) const;
    bool operator<(const NumericType &v) const;

    // change the matrix id of a Numeric Type
    NumericType change_matrix_id(size_t i);

    // change the data id of a Numeric Type
    NumericType change_data_id(size_t i);

    // clear the pool
    static void clear_pool();

    // For checking wheather there are repeated nodes, this will created repeated nodes type if repeatition is detected
    // It will not mark all of the nodes, instead, it will mark the most top level node
    // that is 1. repeated multiple times
    // 2. can conver all of the following repeated nodes
    // it will return true if there are repeated nodes
    // otherwise return false
    void MarkRepeatedNodes(std::map<size_t, size_t> &chosen_repeated_node_map, const bool top_level);

private:
    // make it static
    static std::map<NodeType, char> node_to_op;
};

} // namespace ie
