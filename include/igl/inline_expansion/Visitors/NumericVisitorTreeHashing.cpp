#include "NumericVisitorTreeHashing.hpp"
#include "TreeToFileVisitor.hpp"
#include "TreeToIndexVisitor.hpp"
#include <iostream>
#include <string>
using namespace std;

namespace ie
{
mutex NumericVisitorTreeHashing::can_write_const_map;
mutex NumericVisitorTreeHashing::can_write_operation_map;

NumericVisitorTreeHashing::NumericVisitorTreeHashing() {}

NumericVisitorTreeHashing::NumericVisitorTreeHashing(unsigned int data_size)
{
    this->data_array_operation_ids.resize(data_size);
    this->data_array_used_data_ids.resize(data_size);
    this->data_array_repeated_node_data_ids.resize(data_size);
}

size_t NumericVisitorTreeHashing::write_to_constant_map(const double &const_value)
{
    while (!this->can_write_const_map.try_lock())
    {
        // while some other visit function is writing to constant map
        // we wait
    }
    // when loop exited we have the access to write

    // we need to check if other visitor actually inserted this value mean while
    if (this->seen_consts_position.find(const_value) != this->seen_consts_position.end())
    {
        this->can_write_const_map.unlock();
        return seen_consts_position[const_value];
    }
    // insert it into the map, now we have seen this constant
    this->seen_consts_position.insert({const_value, this->seen_consts.size()});
    this->seen_consts.push_back(const_value);
    size_t position = this->seen_consts.size() - 1;
    // operation done, free the lock
    this->can_write_const_map.unlock();
    // return the correct position in seen_consts
    // this will serve as the left operation for a constant node
    return position;
}

size_t NumericVisitorTreeHashing::write_to_operation_map(const size_t &left_operation_id, const size_t &right_operation_id, const char &current_operation)
{
    while (!this->can_write_operation_map.try_lock())
    {
        // while some other visit function is writing to operation map
        // we wait
    }
    // now claim the operation map write access

    // check if other visotor actually inserted this operation
    if (this->operation_to_id_map.find({{left_operation_id, right_operation_id}, current_operation}) != this->operation_to_id_map.end())
    {
        this->can_write_operation_map.unlock();
        return this->operation_to_id_map[{{left_operation_id, right_operation_id}, current_operation}];
    }
    // cout<<"creating operation "<<this->id_to_operation_map.size()<<": "<<left_operation_id<<" "<<right_operation_id<<" "<<current_operation<<"\n";
    // we initialize this into operation_to_id_map
    this->operation_to_id_map.insert({{{left_operation_id, right_operation_id}, current_operation}, this->id_to_operation_map.size()});
    // push that to id_to_operation_map
    this->id_to_operation_map.push_back({{left_operation_id, right_operation_id}, current_operation});
    // release the lock
    // cout<<"Exiting creation of operation "<<this->id_to_operation_map.size()-1<<"\n";
    size_t position = this->id_to_operation_map.size() - 1;
    this->can_write_operation_map.unlock();
    // return the correct operation id
    return position;
}

void NumericVisitorTreeHashing::visit(NumericType &n, size_t data_position, bool top_level, bool store_position, const map<size_t, size_t> &chosen_repeated_node_map)
{
    if (top_level)
        this->data_array_repeated_node_data_ids[data_position].resize(chosen_repeated_node_map.size()); // we have that many repeated nodes, for each repeated node, it will have its own vector of data ids
    const size_t before_insertion_size = this->data_array_used_data_ids[data_position].size();          // record before inserting ids for this tree, how many data ids we have, this is for the case where it is a repeated tree, then we can just extract those ids
    switch (n.operation)
    {
    case NumericType::Constant:
    {
        if (NumericVisitorTreeHashing::seen_consts_position.find(n.const_value) != NumericVisitorTreeHashing::seen_consts_position.end())
        {
            // if we have seen this const it means it has been recorded not only in
            // seen consts positions, seen consts
            // but also in operations, so we just need to extract that info from operation_to_id_map and store it in
            // data_array_operation_ids
            this->data_array_operation_ids[data_position] = this->operation_to_id_map[{{this->seen_consts_position[n.const_value], 0}, 'c'}];
        }
        else
        {
            // if we have not seen this const
            // we need to initialize everything and create a new operation
            // first we create the position mapping in seen_consts_position
            size_t new_const_position = write_to_constant_map(n.const_value);
            // now write to opeation map
            size_t new_operation_id = write_to_operation_map(new_const_position, 0, 'c');
            // now we record the id of the tree in data_array_operation_ids
            this->data_array_operation_ids[data_position] = new_operation_id;
        }
        break;
    }
    case NumericType::Leaf:
    {
        if (this->operation_to_id_map.find({{n.matrix_id, 0}, 'i'}) != this->operation_to_id_map.end())
        {
            // if we have seen this matrix's stuffs, we need to set the operation, also push the data_id into the array
            this->data_array_operation_ids[data_position] = this->operation_to_id_map[{{n.matrix_id, 0}, 'i'}];
        }
        else
        {
            // create the operation
            size_t new_operation_id = write_to_operation_map(n.matrix_id, 0, 'i');
            // second, set that to the operation id for this position
            this->data_array_operation_ids[data_position] = new_operation_id;
        }
        // now push the data_id
        this->data_array_used_data_ids[data_position].push_back(n.data_id);
        break;
    }
    case NumericType::Add:
    case NumericType::Subtract:
    case NumericType::Multiply:
    case NumericType::Divide:
    {
        // we need to first visit left then visit right
        // this is because we will always store the current operation id in data_array_operation_ids[data_position]
        // we will recursively access it and rewrite it
        // so the order will be: accessing left, write the operation id, parent will access and store that,
        // then access the right, write the operation id, parent access that
        // check if this {{left_operation_id, left_operation_id}, current_operation} exists
        // if not, create and store. If yes, access the operation id and store
        (*this).visit((*NumericType::pool).tree_node_pool[n.left_index], data_position, false, store_position, chosen_repeated_node_map);
        size_t left_operation_id = this->data_array_operation_ids[data_position];
        (*this).visit((*NumericType::pool).tree_node_pool[n.right_index], data_position, false, store_position, chosen_repeated_node_map);
        size_t right_operation_id = this->data_array_operation_ids[data_position];
        // now check if the operation exists
        if (this->operation_to_id_map.find({{left_operation_id, right_operation_id}, n.char_for_operation()}) != this->operation_to_id_map.end())
        {
            // found this operation, store it in data_array_operation_ids
            this->data_array_operation_ids[data_position] = this->operation_to_id_map[{{left_operation_id, right_operation_id}, n.char_for_operation()}];
        }
        else
        {
            // never seen this kind of tree before, add to map
            size_t new_operation_id = write_to_operation_map(left_operation_id, right_operation_id, n.char_for_operation());
            // set the operation in data_array_operation_ids
            this->data_array_operation_ids[data_position] = new_operation_id;
        }
        break;
    }
    case NumericType::Sqrt:
    {
        (*this).visit((*NumericType::pool).tree_node_pool[n.left_index], data_position, false, store_position, chosen_repeated_node_map);
        size_t left_operation_id = this->data_array_operation_ids[data_position];
        // now check if the operation exists
        if (this->operation_to_id_map.find({{left_operation_id, 0}, n.char_for_operation()}) != this->operation_to_id_map.end())
        {
            // found this operation, store it in data_array_operation_ids
            this->data_array_operation_ids[data_position] = this->operation_to_id_map[{{left_operation_id, 0}, n.char_for_operation()}];
        }
        else
        {
            // never seen this kind of tree before, add to map
            size_t new_operation_id = write_to_operation_map(left_operation_id, 0, n.char_for_operation());
            // set the operation in data_array_operation_ids
            this->data_array_operation_ids[data_position] = new_operation_id;
        }
        break;
    }
    }
    // if (top_level){
    //     printf("Top level, operation id is %d, number of repeated entries are %d\n", this->data_array_operation_ids[data_position], chosen_repeated_node_map.size());
    // }
    // if n itself is a repeated node
    if (chosen_repeated_node_map.find(n.self_index) != chosen_repeated_node_map.end())
    {
        const size_t repeated_node_id = chosen_repeated_node_map.at(n.self_index);                // when do we see this repeated node
        const size_t after_insertion_size = this->data_array_used_data_ids[data_position].size(); // how many ids we have now for this tree
        if (this->data_array_repeated_node_data_ids[data_position][repeated_node_id].size() == 0)
        {                                                                                                                                   // if we have not inserted elements before
            this->data_array_repeated_node_data_ids[data_position][repeated_node_id].reserve(after_insertion_size - before_insertion_size); // reserve that many ids
            // std::cout << "Data position " << data_position << " node id " << n.self_index << " repeated id " << repeated_node_id << " used ids: ";
            for (size_t i = before_insertion_size; i < after_insertion_size; i++)
            {
                this->data_array_repeated_node_data_ids[data_position][repeated_node_id].push_back(this->data_array_used_data_ids[data_position][i]); // store the ids into the repeated ids
                // cout << this->data_array_used_data_ids[data_position][i] << " ";
            }
            // std::cout << "\n";
        }
        else
        {
            // std::cout << "Data position " << data_position << " node id " << n.self_index << " repeated id " << repeated_node_id << " used ids: ";
            // for (size_t i = before_insertion_size; i < after_insertion_size; i++)
            // {
            //     cout << this->data_array_used_data_ids[data_position][i] << " ";
            // }
            // std::cout << "\n";
        }
        this->data_array_used_data_ids[data_position].resize(before_insertion_size); // remove those ids
        // now work on inserting this as a repeated node
        if (this->operation_to_id_map.find({{this->data_array_operation_ids[data_position], repeated_node_id}, 'r'}) != this->operation_to_id_map.end())
        {
            this->data_array_operation_ids[data_position] = this->operation_to_id_map[{{this->data_array_operation_ids[data_position], repeated_node_id}, 'r'}];
        }
        else
        {
            size_t new_operation_id = write_to_operation_map(this->data_array_operation_ids[data_position], repeated_node_id, 'r');
            this->data_array_operation_ids[data_position] = new_operation_id;
        }
    }
    if (top_level)
    {
        for (int i = 0; i < chosen_repeated_node_map.size(); i++)
        {
            for (size_t element : this->data_array_repeated_node_data_ids[data_position][i])
                this->data_array_used_data_ids[data_position].push_back(element); // push the repeated data ids at the end
        }
        this->data_array_repeated_node_data_ids[data_position].resize(0);
    }
}

void NumericVisitorTreeHashing::accept(TreeToFileVisitor &visitor)
{
    visitor.visit(*this);
}

void NumericVisitorTreeHashing::accept(TreeToIndexVisitor &visitor)
{
    visitor.visit(*this);
}

} // namespace ie
