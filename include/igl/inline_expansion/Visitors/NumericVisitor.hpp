#pragma once
#include <stddef.h>

namespace ie
{
class NumericType;

class NumericVisitor
{
public:
    virtual void visit(NumericType &n, size_t data_position, bool top_level, bool store_position) = 0;
    NumericVisitor();
    virtual ~NumericVisitor() = 0;
};
} // namespace ie