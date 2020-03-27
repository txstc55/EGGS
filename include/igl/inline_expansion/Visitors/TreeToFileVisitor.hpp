#pragma once
#include <stddef.h>
#include "NumericVisitorTreeHashing.hpp"
#include <inja/inja.hpp>
#include "../md5.h"
#include <iostream>
#include <dlfcn.h>

namespace ie
{
class TreeToFileVisitor
{
private:
public:
    virtual void visit(NumericVisitorTreeHashing &trees) = 0;
    TreeToFileVisitor();
    virtual ~TreeToFileVisitor() = 0;

    // the compiled file name
    // it is a hash of all the hashed function strings, so it is unique I guess
    std::string file_name = "";

    // compile the generated file
    void compile_file();

    virtual void link_functions() = 0;

    // basic function string, it needs to be modified later on to actually fit into function
    std::vector<std::string> operation_strings = {};
    // what matrix id are used here
    std::vector<std::vector<size_t>> operation_matrix_ids = {};

    // operation to string
    std::map<char, std::string> operation_char_to_string = {{'*', "_mm_mul_pd"}, {'+', "_mm_add_pd"}, {'-', "_mm_sub_pd"}, {'/', "_mm_div_pd"}, {'s', "_mm_sqrt_pd"}};

    virtual std::string constant_to_variable_conversion(std::string constant_string);

    // generate all the string in the format of of +(*($, $), 3.1415926) this kind, but needs to be transormed
    // into the form that can be put into a c++ file
    virtual void generate_all_operation_strings(NumericVisitorTreeHashing &trees);

    // convert the previous string to the format of the string that can be used in the file
    virtual std::string generate_usable_string(size_t ind);

    // data to pass into inja for templating
    inja::json data;
};
} // namespace ie