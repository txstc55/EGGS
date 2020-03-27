#include "test_record_util.h"
#include <iostream>
namespace igl
{
IGL_INLINE void write_to_file(std::ofstream &file, std::string message, double time, bool first_call)
{
    if (!first_call)
    {
        file << message << ": " << time << "\n";
        std::cout << message << ": " << time << "\n";
    }
}
} // namespace igl