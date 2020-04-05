#include "TreeToFileVisitorGroupByFunction.hpp"
#include <iostream>
#include <dlfcn.h>

using namespace inja;
using namespace std;

namespace ie
{

TreeToFileVisitorGroupByFunction::TreeToFileVisitorGroupByFunction()
{
    this->data["function_hash_name"] = "";
    this->data["loop_start"] = {};
    this->data["loop_gap"] = {};
    this->data["function_strings"] = {};
    this->data["result_position_start"] = {};
    this->data["loop_length"] = {};
    this->data["matrix_ids"] = {};
    this->data["max_vids"] = 1;
}

TreeToFileVisitorGroupByFunction::TreeToFileVisitorGroupByFunction(size_t gap)
{
    this->data["function_hash_name"] = "";
    this->data["loop_start"] = {};
    this->data["loop_gap"] = {};
    this->data["function_strings"] = {};
    this->data["result_position_start"] = {};
    this->data["loop_length"] = {};
    this->data["matrix_ids"] = {};
    this->data["max_vids"] = 1;
    if (gap > 0)
    {
        this->gap = gap;
    }
}

void TreeToFileVisitorGroupByFunction::visit(NumericVisitorTreeHashing &trees)
{
    this->operation_strings.reserve(trees.id_to_operation_map.size());
    this->operation_matrix_ids.reserve(trees.id_to_operation_map.size());
    generate_all_operation_strings(trees);

    vector<vector<vector<size_t>>> grouped_data_ids;                  // record for each type of tree, what data id are used
    vector<vector<size_t>> grouped_result_position;                   // record for each type of tree, where the result resides
    grouped_data_ids.resize(trees.id_to_operation_map.size());        // there are only that many trees
    grouped_result_position.resize(trees.id_to_operation_map.size()); // there are only that many trees

    size_t loop_index = 0;
    size_t result_index = 0;

    size_t g = this->gap;
    if (g == 0)
        g = trees.data_array_operation_ids.size();
    for (size_t gap_index = 0; gap_index <= trees.data_array_operation_ids.size() / g; gap_index++)
    {

        for (size_t i = gap_index * g; i < ((gap_index + 1) * g) && i < trees.data_array_operation_ids.size(); i++)
        {
            grouped_result_position[trees.data_array_operation_ids[i]].push_back(i);                          // for this type of tree, what position it is in the final output, push it into the group
            grouped_data_ids[trees.data_array_operation_ids[i]].push_back(trees.data_array_used_data_ids[i]); // for this type of tree, what data ids it used, push this vector into the group
        }

        for (size_t i = 0; i < grouped_data_ids.size(); i++)
        {
            if (grouped_result_position[i].size() == 0)
            { // we do nothing
            }
            else
            {
                // setting element for loops
                this->data["loop_start"].push_back(loop_index);                         // where does this loop start
                const size_t operation_ids_size = this->operation_matrix_ids[i].size(); // the number of matrix ids we will need to access

                this->data["loop_gap"].push_back(operation_ids_size); // because there are that many instances for each loop
                // we need to know the max vids used to avoid repeated allocation, this is only for sequential code
                // because we can just allocate that many __m128d at the beginning
                if (this->operation_matrix_ids[i].size() + 1 > this->data["max_vids"])
                {
                    this->data["max_vids"] = this->operation_matrix_ids[i].size() + 1;
                }
                // what matrix id we use
                this->data["matrix_ids"].push_back(this->operation_matrix_ids[i]); // what actual matrix we use for this type of tree
                this->data["result_position_start"].push_back(result_index);       // where will we start store the result for this loop

                // if there are odd number of instances, push back the last one just because we load 2 entries at a time
                if (grouped_data_ids[i].size() % 2 == 1)
                {
                    grouped_data_ids[i].push_back(grouped_data_ids[i][grouped_data_ids[i].size() - 1]);
                    grouped_result_position[i].push_back(grouped_result_position[i][grouped_result_position[i].size() - 1]);
                }

                this->data["loop_length"].push_back(grouped_data_ids[i].size()); // how large is this loop
                result_index += grouped_data_ids[i].size();                      // how many results we stored in this loop
                // where the loop ends
                loop_index += grouped_data_ids[i].size() * operation_ids_size; // how many indices we accessed, this will be used to record where the next loop will start

                // function string
                string usable_string = generate_usable_string(i);        // generate the function string
                this->data["function_strings"].push_back(usable_string); // because we record all the strings
                this->file_name = md5(this->file_name + usable_string);  // the file is only hashed by the function strings not the loop indices, which is bad but I want later on to remove the use of indices so that same structure can have the same file

                // now for our own storage
                this->reordered_data_id.reserve(this->reordered_data_id.size() + operation_ids_size * grouped_data_ids[i].size());   // record in our struct, what the actual data ids we need
                this->reordered_result_position.reserve(this->reordered_result_position.size() + grouped_result_position[i].size()); // record in our struct, what the actual result positions we need
                for (size_t j = 0; j < grouped_data_ids[i].size(); j++)
                    reordered_data_id.insert(reordered_data_id.end(), grouped_data_ids[i][j].begin(), grouped_data_ids[i][j].end());
                reordered_result_position.insert(reordered_result_position.end(), grouped_result_position[i].begin(), grouped_result_position[i].end());
            }
        }
        grouped_result_position.clear();
        grouped_result_position.resize(trees.id_to_operation_map.size());
        grouped_data_ids.clear();
        grouped_data_ids.resize(trees.id_to_operation_map.size());
    }

    this->data["function_hash_name"] = this->file_name;
    Environment env;
    Template function_template = env.parse_template("../include/igl/inline_expansion/GroupByFunction/GroupByFunctionTemplate.txt");
    string result = env.render(function_template, this->data);

    this->file_name = "scramble_grouped_" + this->file_name;
    this->function_hashed_name = this->file_name;
    ofstream foutcpp;
    this->file_name = IE_BIN_DIR + string("/") + this->file_name;
    foutcpp = ofstream(this->file_name + ".cpp");
    foutcpp << result;
    foutcpp.close();
    cout << "Remapped index size: " << this->reordered_data_id.size() << "\n";
    // ofstream map_file;
    // map_file.open("scramble_index.txt");
    // map_file<<reordered_data_id.size()<<" "<<reordered_result_position.size()<<"\n";
    // for (unsigned int i=0; i<reordered_data_id.size(); i++){
    //     map_file<<reordered_data_id[i]<<"\n";
    // }
    // for (unsigned int i=0; i<reordered_result_position.size(); i++){
    //     map_file<<reordered_result_position[i]<<"\n";
    // }
    // map_file.close();
}

void TreeToFileVisitorGroupByFunction::link_functions()
{
#ifdef _WIN32
#define EXTENSION "dll"
#elif defined __unix__
#define EXTENSION "so"
#elif defined __APPLE__
#define EXTENSION "dylib"
#endif
    string lib = this->file_name + string(".") + string(EXTENSION);
    void *handle = dlopen(lib.c_str(), RTLD_NOW);
    if (handle == NULL)
    {
        fprintf(stderr, "dlopen failed: %s\n", dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        *(void **)(&this->Execute_single) = dlsym(handle, (this->function_hashed_name + "_single").c_str());
        if (this->Execute_single == NULL)
        {
            fprintf(stderr, "dlsym for single function failure: %s\n", dlerror());
        }
        *(void **)(&this->Execute_multi) = dlsym(handle, (this->function_hashed_name + "_multi").c_str());
        if (this->Execute_multi == NULL)
        {
            fprintf(stderr, "dlsym for multi function failure: %s\n", dlerror());
        }
    }
}

} // namespace ie