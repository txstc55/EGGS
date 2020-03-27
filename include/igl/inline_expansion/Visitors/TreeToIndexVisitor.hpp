#pragma once
#include <stddef.h>
#include "NumericVisitorTreeHashing.hpp"

namespace ie {

class TreeToIndexVisitor {
public:
    virtual void visit(NumericVisitorTreeHashing& trees) = 0;
    TreeToIndexVisitor();
    virtual ~TreeToIndexVisitor() = 0;
};
}