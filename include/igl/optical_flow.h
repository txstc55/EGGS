#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "igl_inline.h"
#include <string>
#include "mkl_symmetric_solver.h"
#include "inline_expansion/NumericExecutor.hpp"
namespace igl
{
struct OPTICALData
{
    Eigen::MatrixXd image1; // both images are black and white
    Eigen::MatrixXd image2;
    double alpha = 1.0; // weight used in computation
    int method = 0;

    int width, height;

    // all the matrices will be in column major, and turing them into vector or diagonal matrix is quite easy
    Eigen::MatrixXd u; // flow in x direction
    Eigen::MatrixXd v; // flow in y direction
    Eigen::MatrixXd Ix;
    Eigen::MatrixXd Iy;
    Eigen::MatrixXd It;
    Eigen::MatrixXd ubar;
    Eigen::MatrixXd vbar;

    Eigen::SparseMatrix<double> lhs;
    Eigen::VectorXd rhs;

    Eigen::MatrixXd outimage;

    bool first_called = true;

    MKLPardisoData pardiso_data; // for mkl solver

    // for numeric solver
    ie::NumericExecutor ex;                 // the holy executor, our lord and savior
    std::vector<double> result_vector;      // store the result, this is essentially the valuePtr
    std::vector<int> L_outer;               // the outerindexptr
    std::vector<int> L_inner;               // the innerindexptr
    std::vector<std::vector<double>> datas; // datas for executor
    std::vector<double> rhs_vector;         // because we can alsu compute rhs by our way
    ie::NumericExecutor ex_rhs;             // the holy executor, our lord and savior
};

IGL_INLINE void load_image1(OPTICALData &o, std::string file_name);     // load image1
IGL_INLINE void load_image2(OPTICALData &o, std::string file_name);     // load image2, which is the target image
IGL_INLINE void write_image(Eigen::MatrixXd im, std::string file_name); // write out the image, only supports black and white for now

// For those computation, since in the paper it said they should refer to the same pixel, so when accessing the second image, I think it should be image1 plus the flow to get the estimate position of the pixel
IGL_INLINE void compute_ix(OPTICALData &o);   // compute Ex in the paper
IGL_INLINE void compute_iy(OPTICALData &o);   // compute Ey
IGL_INLINE void compute_it(OPTICALData &o);   // compute Et
IGL_INLINE void compute_ubar(OPTICALData &o); // compute the u bar
IGL_INLINE void compute_vbar(OPTICALData &o); // compute the v bar
IGL_INLINE void build_lhs(OPTICALData &o);    // build the lhs of the equation;
IGL_INLINE void build_rhs(OPTICALData &o);    // build the rhs of the equation;
IGL_INLINE void solve_flow(OPTICALData &o);   // solve for u and v

} // namespace igl