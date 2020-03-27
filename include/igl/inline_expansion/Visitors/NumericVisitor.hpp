#pragma once
#include <stddef.h>

namespace ie {
class NumericType;

class NumericVisitor {
public:
    virtual void visit(NumericType& n, size_t data_position) = 0;
    NumericVisitor();
    virtual ~NumericVisitor() = 0;
};
}