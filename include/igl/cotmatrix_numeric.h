#include "igl_inline.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "inline_expansion/NumericExecutor.hpp"

namespace igl
{
struct COTMATRIXData
{                      // Input
    Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
    Eigen::MatrixXi F; // #F by 3/3 list of mesh faces (triangles/tets)

    ie::NumericExecutor ex;                 // the holy executor, our lord and savior
    ie::NumericExecutor ex_l;               // the executor for intermediate resutlts
    std::vector<std::vector<double>> datas; // data data for execution
    bool first_called = true;               // wheather this is the first time we use this struct
    std::vector<int> outerIndex;            // outerindex for the final cotangent matrix
    std::vector<int> innerIndex;            // innerindex for the final cotangent matrix
    std::vector<double> result;             // values for the final cotangent matrix

    std::vector<std::vector<double>> intermediate_datas; // for computing intermediate results
};

IGL_INLINE void squared_edge_lengths_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &V_numeric, const Eigen::MatrixXi &F, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &L); // compute the squared edge length
IGL_INLINE void doublearea_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &l, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &A);                                             // compute the double area
IGL_INLINE void cotmatrix_entries_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &V_numeric, const Eigen::MatrixXi &F, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &C);    // compute the entries of the cotangent matrix
IGL_INLINE void cotmatrix_numeric(COTMATRIXData &cd);                                                                                                                                                                 // putting them into the sparse matrix

IGL_INLINE void cotmatrix_entries_numeric_intermediate(COTMATRIXData &cd, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &C); // compute the entries of the cotangent matrix but we will compute some data intermediately not in one shot
IGL_INLINE void cotmatrix_numeric_intermediate(COTMATRIXData &cd);                                                                            // putting them into the sparse matrix

} // namespace igl