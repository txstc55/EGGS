#include "cotmatrix_numeric.h"
#include "inline_expansion/utils.hpp"
namespace igl
{

IGL_INLINE void squared_edge_lengths_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &V_numeric, const Eigen::MatrixXi &F, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &L)
{
    using namespace std;
    const int m = F.rows();
    L.resize(m, 3);
    for (int i = 0; i < m; i++)
    {
        L(i, 0) = (V_numeric(F(i, 1), 0) - V_numeric(F(i, 2), 0)) * (V_numeric(F(i, 1), 0) - V_numeric(F(i, 2), 0)) + (V_numeric(F(i, 1), 1) - V_numeric(F(i, 2), 1)) * (V_numeric(F(i, 1), 1) - V_numeric(F(i, 2), 1)) + (V_numeric(F(i, 1), 2) - V_numeric(F(i, 2), 2)) * (V_numeric(F(i, 1), 2) - V_numeric(F(i, 2), 2));
        L(i, 1) = (V_numeric(F(i, 2), 0) - V_numeric(F(i, 0), 0)) * (V_numeric(F(i, 2), 0) - V_numeric(F(i, 0), 0)) + (V_numeric(F(i, 2), 1) - V_numeric(F(i, 0), 1)) * (V_numeric(F(i, 2), 1) - V_numeric(F(i, 0), 1)) + (V_numeric(F(i, 2), 2) - V_numeric(F(i, 0), 2)) * (V_numeric(F(i, 2), 2) - V_numeric(F(i, 0), 2));
        L(i, 2) = (V_numeric(F(i, 0), 0) - V_numeric(F(i, 1), 0)) * (V_numeric(F(i, 0), 0) - V_numeric(F(i, 1), 0)) + (V_numeric(F(i, 0), 1) - V_numeric(F(i, 1), 1)) * (V_numeric(F(i, 0), 1) - V_numeric(F(i, 1), 1)) + (V_numeric(F(i, 0), 2) - V_numeric(F(i, 1), 2)) * (V_numeric(F(i, 0), 2) - V_numeric(F(i, 1), 2));
    }
}

IGL_INLINE void doublearea_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &l, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &A)
{
    A.resize(l.rows(), 1);
    for (int i = 0; i < l.rows(); i++)
    {
        A(i, 0) = 0.5 * ((l(i, 0) + (l(i, 1) + l(i, 2))) *
                         (l(i, 2) - (l(i, 0) - l(i, 1))) *
                         (l(i, 2) + (l(i, 0) - l(i, 1))) *
                         (l(i, 0) + (l(i, 1) - l(i, 2))))
                            .sqrt();
    }
}

IGL_INLINE void cotmatrix_entries_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &V_numeric, const Eigen::MatrixXi &F, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &C)
{
    using namespace std;
    using namespace Eigen;
    Matrix<ie::NumericType, Dynamic, Dynamic> l2;
    // Compute the squared edge length
    squared_edge_lengths_numeric(V_numeric, F, l2);
    // Compute Edge lengths
    Matrix<ie::NumericType, Dynamic, Dynamic> l;
    l.resize(l2.rows(), 3);
    for (int i = 0; i < l.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
            l(i, j) = l2(i, j).sqrt();
    }

    Matrix<ie::NumericType, Dynamic, Dynamic> dblA;
    doublearea_numeric(l, dblA);
    C.resize(F.rows(), 3);
    for (int i = 0; i < F.rows(); i++)
    {
        C(i, 0) = (l2(i, 1) + l2(i, 2) - l2(i, 0)) / dblA(i) / 4.0;
        C(i, 1) = (l2(i, 2) + l2(i, 0) - l2(i, 1)) / dblA(i) / 4.0;
        C(i, 2) = (l2(i, 0) + l2(i, 1) - l2(i, 2)) / dblA(i) / 4.0;
    }
}

IGL_INLINE void cotmatrix_entries_numeric_intermediate(COTMATRIXData &cd, Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &C)
{
    using namespace std;
    using namespace Eigen;
    cd.datas.resize(1);
    cd.intermediate_datas.resize(1);
    cd.intermediate_datas[0].resize(cd.V.rows() * cd.V.cols());
    Matrix<ie::NumericType, Dynamic, Dynamic> l2;
    Matrix<ie::NumericType, Dynamic, Dynamic> V_numeric;
    V_numeric.resize(cd.V.rows(), cd.V.cols());
    V_numeric = ie::to_dense_numeric(cd.V, 0);
    // Compute the squared edge length
    squared_edge_lengths_numeric(V_numeric, cd.F, l2);
    cd.ex_l = ie::NumericExecutor(l2, 0);      // from this point on, everything else depends on the value of l2, ex_l takes the value in V as input
    cd.datas[0].resize(l2.rows() * l2.cols()); // we will use l2 as the input for the next executor, so we can just put the results in data

    MatrixXd tmp_l2; // now make l2 the input instead of V
    tmp_l2.resize(l2.rows(), l2.cols());
    Matrix<ie::NumericType, Dynamic, Dynamic> tmp_l2_numeric; // routines, probably should do the resize in utils instead
    tmp_l2_numeric.resize(l2.rows(), l2.cols());
    tmp_l2_numeric = ie::to_dense_numeric(tmp_l2, 0);

    // Compute Edge lengths
    Matrix<ie::NumericType, Dynamic, Dynamic> l;
    l.resize(l2.rows(), 3);
    for (int i = 0; i < l.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
            l(i, j) = tmp_l2_numeric(i, j).sqrt();
    }

    Matrix<ie::NumericType, Dynamic, Dynamic> dblA;
    doublearea_numeric(l, dblA);
    C.resize(cd.F.rows(), 3);
    for (int i = 0; i < cd.F.rows(); i++)
    {
        C(i, 0) = (tmp_l2_numeric(i, 1) + tmp_l2_numeric(i, 2) - tmp_l2_numeric(i, 0)) / dblA(i) / 4.0;
        C(i, 1) = (tmp_l2_numeric(i, 2) + tmp_l2_numeric(i, 0) - tmp_l2_numeric(i, 1)) / dblA(i) / 4.0;
        C(i, 2) = (tmp_l2_numeric(i, 0) + tmp_l2_numeric(i, 1) - tmp_l2_numeric(i, 2)) / dblA(i) / 4.0;
    }
}

IGL_INLINE void cotmatrix_numeric(COTMATRIXData &cd)
{
    using namespace Eigen;
    using namespace std;
    if (cd.first_called)
    {
        cd.datas.resize(1);
        cd.datas[0].resize(cd.V.rows() * cd.V.cols());
        SparseMatrix<ie::NumericType> L;
        L.resize(cd.V.rows(), cd.V.rows());
        Matrix<int, Dynamic, 2> edges;
        L.reserve(10 * cd.V.rows());
        edges.resize(3, 2);
        edges << 1, 2,
            2, 0,
            0, 1;
        Matrix<ie::NumericType, Dynamic, Dynamic> C;
        Matrix<ie::NumericType, Dynamic, Dynamic> V_numeric;
        V_numeric.resize(cd.V.rows(), cd.V.cols());
        V_numeric = ie::to_dense_numeric(cd.V, 0);
        cotmatrix_entries_numeric(V_numeric, cd.F, C);

        vector<Triplet<ie::NumericType>> IJV;
        IJV.reserve(cd.F.rows() * edges.rows() * 4);
        // Loop over triangles
        for (int i = 0; i < cd.F.rows(); i++)
        {
            // loop over edges of element
            for (int e = 0; e < edges.rows(); e++)
            {
                int source = cd.F(i, edges(e, 0));
                int dest = cd.F(i, edges(e, 1));
                IJV.push_back(Triplet<ie::NumericType>(source, dest, C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(dest, source, C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(source, source, ie::NumericType(0.0) - C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(dest, dest, ie::NumericType(0.0) - C(i, e)));
            }
        }
        L.setFromTriplets(IJV.begin(), IJV.end());
        L.makeCompressed();
        cd.ex = ie::NumericExecutor(L, 0);
        cd.result.resize(L.nonZeros());
        cd.first_called = false;
    }
    cd.datas[0].assign(cd.V.data(), cd.V.data() + cd.V.rows() * cd.V.cols());
    cd.ex.ExecuteMulti(cd.datas, cd.result);
}

IGL_INLINE void cotmatrix_numeric_intermediate(COTMATRIXData &cd)
{
    using namespace Eigen;
    using namespace std;
    if (cd.first_called)
    {
        SparseMatrix<ie::NumericType> L;
        L.resize(cd.V.rows(), cd.V.rows());
        Matrix<int, Dynamic, 2> edges;
        L.reserve(10 * cd.V.rows());
        edges.resize(3, 2);
        edges << 1, 2,
            2, 0,
            0, 1;
        Matrix<ie::NumericType, Dynamic, Dynamic> C;
        cotmatrix_entries_numeric_intermediate(cd, C);

        vector<Triplet<ie::NumericType>> IJV;
        IJV.reserve(cd.F.rows() * edges.rows() * 4);
        // Loop over triangles
        for (int i = 0; i < cd.F.rows(); i++)
        {
            // loop over edges of element
            for (int e = 0; e < edges.rows(); e++)
            {
                int source = cd.F(i, edges(e, 0));
                int dest = cd.F(i, edges(e, 1));
                IJV.push_back(Triplet<ie::NumericType>(source, dest, C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(dest, source, C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(source, source, ie::NumericType(0.0) - C(i, e)));
                IJV.push_back(Triplet<ie::NumericType>(dest, dest, ie::NumericType(0.0) - C(i, e)));
            }
        }
        L.setFromTriplets(IJV.begin(), IJV.end());
        L.makeCompressed();
        cd.ex = ie::NumericExecutor(L, 0);
        cd.result.resize(L.nonZeros());
        cd.first_called = false;
    }
    cd.intermediate_datas[0].assign(cd.V.data(), cd.V.data() + cd.V.rows() * cd.V.cols());
    cd.ex_l.ExecuteMulti(cd.intermediate_datas, cd.datas[0]);
    cd.ex.ExecuteMulti(cd.datas, cd.result);
}

} // namespace igl