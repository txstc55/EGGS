#include "TreeToFileVisitor.hpp"
using namespace std;
#include <sys/types.h>
#include <sys/wait.h>
#include <dlfcn.h>
#include <unistd.h>
#include <iostream>

namespace ie
{
TreeToFileVisitor::TreeToFileVisitor() {}
TreeToFileVisitor::~TreeToFileVisitor() {}

string TreeToFileVisitor::constant_to_variable_conversion(string constant_string)
{
    return "__m128d{" + constant_string + ", " + constant_string + "}";
}

void TreeToFileVisitor::generate_all_operation_strings(NumericVisitorTreeHashing &trees)
{
    for (unsigned int i = 0; i < trees.id_to_operation_map.size(); i++)
    {
        switch (trees.id_to_operation_map[i].second)
        {
        case 'c':
        {
            this->operation_strings.push_back(constant_to_variable_conversion(to_string(trees.seen_consts[trees.id_to_operation_map[i].first.first])));
            this->operation_matrix_ids.push_back({});
            break;
        }
        case 'i':
        {
            this->operation_strings.push_back("$");
            this->operation_matrix_ids.push_back({trees.id_to_operation_map[i].first.first});
            break;
        }
        case 's':
        {
            string function_string = "";
            function_string += this->operation_char_to_string[trees.id_to_operation_map[i].second];
            function_string += "(";
            function_string += this->operation_strings[trees.id_to_operation_map[i].first.first] + ")";
            this->operation_strings.push_back(function_string);
            vector<size_t> tmp(this->operation_matrix_ids[trees.id_to_operation_map[i].first.first]);
            this->operation_matrix_ids.push_back(tmp);
            break;
        }
        case 'r':
        {
            string function_string = "";
            function_string += this->operation_char_to_string[trees.id_to_operation_map[i].second];
            function_string += "(";
            function_string += this->operation_strings[trees.id_to_operation_map[i].first.first] + ")";
            this->operation_strings.push_back(function_string);
            vector<size_t> tmp(this->operation_matrix_ids[trees.id_to_operation_map[i].first.first]);
            this->operation_matrix_ids.push_back(tmp);
            break;
        }
        default:
        {
            string function_string = "";
            function_string += this->operation_char_to_string[trees.id_to_operation_map[i].second];
            function_string += "(";
            function_string += this->operation_strings[trees.id_to_operation_map[i].first.first] + ", " + this->operation_strings[trees.id_to_operation_map[i].first.second] + ")";
            this->operation_strings.push_back(function_string);
            vector<size_t> tmp(this->operation_matrix_ids[trees.id_to_operation_map[i].first.first]);
            tmp.insert(tmp.end(), this->operation_matrix_ids[trees.id_to_operation_map[i].first.second].begin(), this->operation_matrix_ids[trees.id_to_operation_map[i].first.second].end());
            this->operation_matrix_ids.push_back(tmp);
            break;
        }
        }
    }
}

string TreeToFileVisitor::generate_usable_string(size_t ind)
{
    string final_output = "";
    unsigned int count = 0;
    for (unsigned int i = 0; i < this->operation_strings[ind].size(); i++)
    {
        if (this->operation_strings[ind][i] == '$')
        {
            final_output += "v" + to_string(count);
            count++;
        }
        else
        {
            final_output += this->operation_strings[ind][i];
        }
    }
    return final_output;
}

void TreeToFileVisitor::compile_file()
{
#ifdef _WIN32
#define EXTENSION "dll"
#elif defined __unix__
#define EXTENSION "so"
#elif defined __APPLE__
#define EXTENSION "dylib"
#endif

    // compile the file
    static const string ie_bin_dir = IE_BIN_DIR;
    static const string ie_cxx_compiler = IE_CXX_COMPILER;
    static const string tbb_dir = IE_TBB_DIR;
    static const string tbb_lib_dir = IE_TBB_LIB;

    char *compiler = new char[string("clang++").length() + 1];
    strcpy(compiler, string("clang++").c_str());
    char *filename = new char[string(this->file_name + ".cpp").length() + 1];
    strcpy(filename, string(this->file_name + ".cpp").c_str());
    char *rpath = new char[(string("-Wl,-rpath,") + IE_TBB_LIB + string("/build")).length() + 1];
    strcpy(rpath, (string("-Wl,-rpath,") + IE_TBB_LIB + string("/build")).c_str());
    char *include = new char[(string("-I") + IE_TBB_DIR + string("/include")).length() + 1];
    strcpy(include, (string("-I") + IE_TBB_DIR + string("/include")).c_str());
    char *link = new char[(string("-L") + IE_TBB_LIB + string("/..")).length() + 1];
    strcpy(link, (string("-L") + IE_TBB_LIB + string("/..")).c_str());
    char *outlib = new char[("-o" + this->file_name + "." + EXTENSION).length() + 1];
    strcpy(outlib, ("-o" + this->file_name + "." + EXTENSION).c_str());
    char *tbb_static = new char[string("-ltbb_static").length() + 1];
    strcpy(tbb_static, string("-ltbb_static").c_str());
    char *avx2 = new char[string("-mavx2").length() + 1];
    strcpy(avx2, string("-mavx2").c_str());
    char *sse4 = new char[7];
    strcpy(sse4, string("-msse4").c_str());
    char *shared = new char[8];
    strcpy(shared, string("-shared").c_str());
    char *fpic = new char[6];
    strcpy(fpic, string("-fPIC").c_str());
    char *std11 = new char[string("-std=c++11").length() + 1];
    strcpy(std11, string("-std=c++11").c_str());
    char *o3 = new char[7];
    strcpy(o3, string("-Ofast").c_str());
    char *native = new char[string("-march=native").length() + 1];
    strcpy(native, string("-march=native").c_str());
    string command = string("clang++") + " " + this->file_name + ".cpp -Wl,-rpath," + IE_TBB_LIB + "/build -I" + IE_TBB_DIR + "/include -L" + IE_TBB_LIB + "/.." + " -ltbb_static -mavx2 -msse4 -o " + this->file_name + "." + EXTENSION + " -shared -fPIC -std=c++11 -Ofast -march=native";
    char *args[] = {compiler, filename, rpath, include, link, tbb_static, avx2, sse4, outlib, shared, fpic, std11, o3, native, NULL};
    cout << "Forking a child process\n";
    pid_t pid;
    pid = fork();
    if (pid == -1)
    {
        cout << "Failed to fork child\n";
    }
    else if (pid == 0)
    {
        cout << command << "\n";
        cout << "Compiling generated file\n";
        execvp(args[0], args);
    }
    else
    {
        wait(NULL);
    }
}

} // namespace ie