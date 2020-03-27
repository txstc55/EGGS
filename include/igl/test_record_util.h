#pragma once
#include "igl_inline.h"
#include <iostream>
#include <fstream>

namespace igl
{
IGL_INLINE void write_to_file(std::ofstream &file, std::string message, double time, bool first_call);
}