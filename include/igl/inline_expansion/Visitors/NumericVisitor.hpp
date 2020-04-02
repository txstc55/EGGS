#pragma once
#include <stddef.h>
#include <set>
namespace ie
{
class NumericType;

class NumericVisitor
{
public:
    virtual void visit(NumericType &n, size_t data_position, const bool top_level, const bool store_position, const std::set<size_t> &chosen_repeated_node) = 0;
    NumericVisitor();
    virtual ~NumericVisitor() = 0;
};
} // namespace ie