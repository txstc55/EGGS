#pragma once
#include <Eigen/Sparse>
#include "igl_inline.h"
#include <fstream>
#include <iostream>
namespace igl
{
IGL_INLINE void get_triplet_value_order_col_major(const std::vector<Eigen::Triplet<double>> &trip, std::vector<int> &order);                    // assuming the order of triplet remains the same, find out where in the colum major matrix's value pointer should we put each triplet's value
IGL_INLINE void get_triplet_value_order_row_major(const std::vector<Eigen::Triplet<double>> &trip, std::vector<int> &order);                    // assuming the order of triplet remains the same, find out where in the colum major matrix's value pointer should we put each triplet's value
IGL_INLINE void put_triplet_value(const std::vector<Eigen::Triplet<double>> &trip, const std::vector<int> &order, std::vector<double> &values); // put the values into the corresponding order
IGL_INLINE void export_mtx(const Eigen::SparseMatrix<double> &m, std::string file_name);                                                        //export to mtx format
} // namespace igl