#include "ie_helper.hpp"

namespace igl
{
IGL_INLINE void get_triplet_value_order_col_major(const std::vector<Eigen::Triplet<double>> &trip, std::vector<int> &order)
{
    std::vector<std::vector<int>> T(trip.size());
    for (unsigned i = 0; i < trip.size(); ++i)
    {
        T[i].resize(3);
        T[i][0] = trip[i].col();
        T[i][1] = trip[i].row();
        T[i][2] = i;
    }

    std::sort(T.begin(), T.end());
    order.resize(T.size());

    for (unsigned int i = 0; i < T.size(); i++)
    {
        order[T[i][2]] = i;
    }
}

IGL_INLINE void get_triplet_value_order_row_major(const std::vector<Eigen::Triplet<double>> &trip, std::vector<int> &order)
{
    std::vector<std::vector<int>> T(trip.size());
    for (unsigned i = 0; i < trip.size(); ++i)
    {
        T[i].resize(3);
        T[i][0] = trip[i].row();
        T[i][1] = trip[i].col();
        T[i][2] = i;
    }

    std::sort(T.begin(), T.end());
    order.resize(T.size());

    for (unsigned int i = 0; i < T.size(); i++)
    {
        order[T[i][2]] = i;
    }
}

IGL_INLINE void put_triplet_value(const std::vector<Eigen::Triplet<double>> &trip, const std::vector<int> &order, std::vector<double> &values)
{
    for (unsigned int i = 0; i < trip.size(); i++)
    {
        values[order[i]] = trip[i].value();
    }
}

IGL_INLINE void export_mtx(const Eigen::SparseMatrix<double> &m, std::string file_name)
{
    std::ofstream matrix_file(file_name);
    matrix_file << "%%MatrixMarket matrix coordinate real general\n";
    matrix_file << m.rows() << " " << m.cols() << " " << m.nonZeros() << "\n";
    for (int k = 0; k < m.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it)
        {
            matrix_file << it.row() + 1 << " " << it.col() + 1 << " " << it.value() << "\n";
        }
    }
    matrix_file.close();
    std::cout << "First time writing file done\n";
}

} // namespace igl