#include "mkl_symmetric_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

namespace igl
{

IGL_INLINE void pardiso_init(MKLPardisoData &m)
{
    for (int i = 0; i < 64; i++)
    {
        m.iparm[i] = 0;
    }
    m.iparm[0] = 1;   /* No solver default */
    m.iparm[1] = 2;   /* Fill-in reordering from METIS */
    m.iparm[3] = 0;   /* No iterative-direct algorithm */
    m.iparm[4] = 0;   /* No user fill-in reducing permutation */
    m.iparm[5] = 0;   /* Write solution into x */
    m.iparm[6] = 0;   /* Not in use */
    m.iparm[7] = 2;   /* Max numbers of iterative refinement steps */
    m.iparm[8] = 0;   /* Not in use */
    m.iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    m.iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    m.iparm[11] = 0;  /* Not in use */
    m.iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try m.iparm[12] = 1 in case of inappropriate accuracy */
    m.iparm[13] = 0;  /* Output: Number of perturbed pivots */
    m.iparm[14] = 0;  /* Not in use */
    m.iparm[15] = 0;  /* Not in use */
    m.iparm[16] = 0;  /* Not in use */
    m.iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    m.iparm[18] = -1; /* Output: Mflops for LU factorization */
    m.iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    m.iparm[34] = 1;  /* zero index based csr format */
    m.maxfct = 1;     /* Maximum number of numerical factorizations. */
    m.mnum = 1;       /* Which factorization to use. */
    m.msglvl = 0;     /* Print statistical information in file 0 for disable, 1 for enable */
    m.error = 0;      /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < 64; i++)
    {
        m.pt[i] = 0;
    }
}

IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, int *outer, int *inner, double *value, int n)
{
    m.n = n;
    m.outerIndexPtr = outer;
    m.innerIndexPtr = inner;
    m.valuePtr = value;
}

IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, Eigen::SparseMatrix<double> &mat)
{
    pardiso_support_matrix(m, mat.outerIndexPtr(), mat.innerIndexPtr(), mat.valuePtr(), mat.rows());
}

IGL_INLINE void pardiso_support_matrix(MKLPardisoData &m, Eigen::SparseMatrix<double, Eigen::RowMajor> &mat)
{
    pardiso_support_matrix(m, mat.outerIndexPtr(), mat.innerIndexPtr(), mat.valuePtr(), mat.rows());
}

IGL_INLINE void pardiso_support_value(MKLPardisoData &m, double *value)
{
    m.valuePtr = value;
}

IGL_INLINE void pardiso_symbolic_factor(MKLPardisoData &m)
{
    m.phase = 11;
    PARDISO(m.pt, &m.maxfct, &m.mnum, &m.mtype, &m.phase,
            &m.n, m.valuePtr, m.outerIndexPtr, m.innerIndexPtr, &m.idum, &m.nrhs, m.iparm, &m.msglvl, &m.ddum, &m.ddum, &m.error);
    if (m.error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", m.error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors = %d", m.iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", m.iparm[18]);
}

IGL_INLINE void pardiso_numeric_factor(MKLPardisoData &m)
{
    m.phase = 22;
    PARDISO(m.pt, &m.maxfct, &m.mnum, &m.mtype, &m.phase,
            &m.n, m.valuePtr, m.outerIndexPtr, m.innerIndexPtr, &m.idum, &m.nrhs, m.iparm, &m.msglvl, &m.ddum, &m.ddum, &m.error);
    if (m.error != 0)
    {
        printf("\nERROR during numerical factorization: %d", m.error);
        exit(2);
    }
    printf("\nFactorization completed ... ");
}

IGL_INLINE void pardiso_solve(MKLPardisoData &m, double *x, double *b)
{
    m.phase = 33;
    PARDISO(m.pt, &m.maxfct, &m.mnum, &m.mtype, &m.phase,
            &m.n, m.valuePtr, m.outerIndexPtr, m.innerIndexPtr, &m.idum, &m.nrhs, m.iparm, &m.msglvl, b, x, &m.error);
    if (m.error != 0)
    {
        printf("\nERROR during solution: %d", m.error);
        exit(3);
    }
}

} // namespace igl