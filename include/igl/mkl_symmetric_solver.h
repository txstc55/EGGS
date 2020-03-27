#pragma once
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

#include "igl_inline.h"

#include <Eigen/Sparse>
namespace igl
{
struct MKLPardisoData
{
    MKL_INT n;                                  // the number of equations
    MKL_INT mtype = -2;                         // Real symmetric matrix
    double res, res0;                           // residuals i guess
    MKL_INT nrhs = 1;                           // Number of right hand sides.
    void *pt[64];                               //solving double
    MKL_INT iparm[64];                          //paramaters used in pardiso
    MKL_INT maxfct, mnum, phase, error, msglvl; // miscs used during the solve
    double ddum;                                /* Double dummy */
    MKL_INT idum;                               /* Integer dummy. */

    // matrix csr format
    int *outerIndexPtr;
    int *innerIndexPtr;
    double *valuePtr;
};

IGL_INLINE void pardiso_init(MKLPardisoData &m);                                                              // init the paramaters for a general matrix
IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, Eigen::SparseMatrix<double> &mat);                  // support the matrix in eigen
IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, Eigen::SparseMatrix<double, Eigen::RowMajor> &mat); // support the matrix in eigen
IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, int *outer, int *inner, double *value, int n);      // support the matrix in csr format
IGL_INLINE void pardiso_support_value(MKLPardisoData &m, double *value);                                      // keeping the structure but support another value array
IGL_INLINE void pardiso_symbolic_factor(MKLPardisoData &m);                                                   // symbolic facatorization
IGL_INLINE void pardiso_numeric_factor(MKLPardisoData &m);                                                    // numerical factorization
IGL_INLINE void pardiso_solve(MKLPardisoData &m, double *x, double *b);                                       // compute the result into x, need the support for residuals
} // namespace igl
