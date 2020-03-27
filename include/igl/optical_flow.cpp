#include "optical_flow.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include "test_record_util.h"
#include "Timer.h"
#include "inline_expansion/utils.hpp"

namespace igl
{

std::ofstream result_file; // the file to write the test result

// convert double to char and char to double for stb image
unsigned char double_2_unsignedchar(const double d)
{
    return round(std::max(std::min(255., d), 0.));
}
double unsignedchar_2_double(const unsigned char c)
{
    return (double)c;
}

// return value if not out of bound
ie::NumericType return_matrix_value_numeric(Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> &m, const int i, const int j)
{
    if (i == -1 || j == -1 || i == m.rows() || j == m.cols())
    {
        return 0;
    }
    else
    {
        return m(i, j);
    }
}
double return_matrix_value(const Eigen::MatrixXd &m, const int i, const int j)
{
    if (i == -1 || j == -1 || i == m.rows() || j == m.cols())
    {
        return 0;
    }
    else
    {
        return m(i, j);
    }
}

// write an eigen matrix to file
IGL_INLINE void write_image(Eigen::MatrixXd im, std::string file_name)
{
    const int w = im.cols();                          // Image width
    const int h = im.rows();                          // Image height
    const int comp = 1;                               // 3 Channels Red, Green, Blue, Alpha
    const int stride_in_bytes = w * comp;             // Lenght of one row in bytes
    std::vector<unsigned char> data(w * h * comp, 0); // The image itself;

    for (unsigned wi = 0; wi < w; ++wi)
    {
        for (unsigned hi = 0; hi < h; ++hi)
        {
            data[(hi * w) + wi] = double_2_unsignedchar(im(hi, wi));
        }
    }
    stbi_write_png(file_name.c_str(), w, h, comp, data.data(), stride_in_bytes);
}

// load an image to optical data. It will always be black and white
IGL_INLINE void load_image1(OPTICALData &o, std::string file_name)
{
    int n;
    unsigned char *data = stbi_load(file_name.c_str(), &o.width, &o.height, &n, 0);
    std::cout << n << "<======== number of components\n";
    o.image1.resize(o.height, o.width);
    for (unsigned wi = 0; wi < o.width; ++wi)
    {
        for (unsigned hi = 0; hi < o.height; ++hi)
        {
            for (int i = 0; i < n; i++)
            {
                o.image1(hi, wi) += unsignedchar_2_double(data[(hi * o.width + wi) * n + i]);
            }
            o.image1(hi, wi) /= (n * 1.0);
        }
    }

    if (o.first_called)
    {
        o.Ix.resize(o.height, o.width);
        o.Iy.resize(o.height, o.width);
        o.It.resize(o.height, o.width);
        o.ubar.resize(o.height, o.width);
        o.vbar.resize(o.height, o.width);
        o.u.resize(o.height, o.width);
        o.v.resize(o.height, o.width);
        o.lhs.resize(o.height * o.width * 2, o.height * o.width * 2);
        for (int i = 0; i < o.height; i++)
        {
            for (int j = 0; j < o.width; j++)
            {
                o.Ix(i, j) = 0.0;
                o.Iy(i, j) = 0.0;
                o.It(i, j) = 0.0;
                o.ubar(i, j) = 0.0;
                o.vbar(i, j) = 0.0;
                o.u(i, j) = 0.0;
                o.v(i, j) = 0.0;
            }
        }
    }
}

// load the destination image
IGL_INLINE void load_image2(OPTICALData &o, std::string file_name)
{
    int n;
    unsigned char *data = stbi_load(file_name.c_str(), &o.width, &o.height, &n, 0);
    o.image2.resize(o.height, o.width);
    for (unsigned wi = 0; wi < o.width; ++wi)
    {
        for (unsigned hi = 0; hi < o.height; ++hi)
        {
            for (int i = 0; i < n; i++)
            {
                o.image2(hi, wi) += unsignedchar_2_double(data[(hi * o.width + wi) * n + i]);
            }
            o.image2(hi, wi) /= (n * 1.0);
        }
    }
}

IGL_INLINE void compute_ix(OPTICALData &o)
{
    for (int w = 0; w < o.width - 1; w++)
    {
        for (int h = 0; h < o.height - 1; h++)
        {
            o.Ix(h, w) = (o.image1(h, w + 1) - o.image1(h, w) + o.image1(h + 1, w + 1) - o.image1(h + 1, w) + o.image2(h, w + 1) - o.image2(h, w) + o.image2(h + 1, w + 1) - o.image2(h + 1, w)) / 4.0;
        }
    }
    // set the boundary now
    for (int w = 0; w < o.width - 1; w++)
    {
        o.Ix(o.height - 1, w) = 0.0;
    }
    for (int h = 0; h < o.height; h++)
    {
        o.Ix(h, o.width - 1) = 0.0;
    }
}

IGL_INLINE void compute_iy(OPTICALData &o)
{
    for (int w = 0; w < o.width - 1; w++)
    {
        for (int h = 0; h < o.height - 1; h++)
        {
            o.Iy(h, w) = (o.image1(h + 1, w) - o.image1(h, w) + o.image1(h + 1, w + 1) - o.image1(h, w + 1) + o.image2(h + 1, w) - o.image2(h, w) + o.image2(h + 1, w + 1) - o.image2(h, w + 1)) / 4.0;
        }
    }
    // set the boundary now
    for (int h = 0; h < o.height - 1; h++)
    {
        o.Iy(h, o.width - 1) = 0.0;
    }

    for (int w = 0; w < o.width; w++)
    {
        o.Iy(o.height - 1, w) = 0.0;
    }
}

IGL_INLINE void compute_it(OPTICALData &o)
{
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            o.It(h, w) = (return_matrix_value(o.image2, h, w) - return_matrix_value(o.image1, h, w) + return_matrix_value(o.image2, h + 1, w) - return_matrix_value(o.image1, h + 1, w) + return_matrix_value(o.image2, h, w + 1) - return_matrix_value(o.image1, h, w + 1) + return_matrix_value(o.image2, h + 1, w + 1) - return_matrix_value(o.image1, h + 1, w + 1)) / 4.0;
        }
    }
}

IGL_INLINE void compute_ubar(OPTICALData &o)
{
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            o.ubar(h, w) = (return_matrix_value(o.u, h + 1, w) + return_matrix_value(o.u, h, w + 1) + return_matrix_value(o.u, h - 1, w) + return_matrix_value(o.u, h, w - 1)) / 6.0 + (return_matrix_value(o.u, h - 1, w - 1) + return_matrix_value(o.u, h - 1, w + 1) + return_matrix_value(o.u, h + 1, w - 1) + return_matrix_value(o.u, h + 1, w + 1)) / 12.0;
        }
    }
}

IGL_INLINE void compute_vbar(OPTICALData &o)
{
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            o.vbar(h, w) = (return_matrix_value(o.v, h + 1, w) + return_matrix_value(o.v, h, w + 1) + return_matrix_value(o.v, h - 1, w) + return_matrix_value(o.v, h, w - 1)) / 6.0 + (return_matrix_value(o.v, h - 1, w - 1) + return_matrix_value(o.v, h - 1, w + 1) + return_matrix_value(o.v, h + 1, w - 1) + return_matrix_value(o.v, h + 1, w + 1)) / 12.0;
        }
    }
}

IGL_INLINE void build_lhs(OPTICALData &o)
{
    result_file << "START BUILDING LEFT HAND SIDE USING EIGEN\n";
    igl::Timer t;
    t.start();
    compute_ix(o);
    compute_iy(o);
    compute_it(o);
    t.stop();
    write_to_file(result_file, "COMPUTE IX, IY, IT", t.getElapsedTimeInMicroSec(), o.first_called);
    t.start();
    double a2 = std::pow(o.alpha, 2);
    Eigen::MatrixXd x_squared = o.Ix.cwiseProduct(o.Ix) + Eigen::MatrixXd::Constant(o.height, o.width, a2);
    Eigen::MatrixXd y_squared = o.Iy.cwiseProduct(o.Iy) + Eigen::MatrixXd::Constant(o.height, o.width, a2);
    Eigen::MatrixXd xy = o.Ix.cwiseProduct(o.Iy);
    t.stop();
    write_to_file(result_file, "EIGEN COMPUTING ENTRIES OF SPARSE MATRIX", t.getElapsedTimeInMicroSec(), o.first_called);
    t.start();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(o.height * o.width * 3);
    int index = 0;
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            trip.push_back(Eigen::Triplet<double>(index, index, x_squared(h, w)));
            index++;
        }
    }
    int index2 = index;
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            trip.push_back(Eigen::Triplet<double>(index, index, y_squared(h, w)));
            index++;
        }
    }
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            trip.push_back(Eigen::Triplet<double>(index2, index2 - o.height * o.width, xy(h, w)));
            index2++;
        }
    }
    o.lhs.setFromTriplets(trip.begin(), trip.end());
    o.lhs.makeCompressed();
    t.stop();
    write_to_file(result_file, "EIGEN ASSEMBLE SPARSE MATRIX", t.getElapsedTimeInMicroSec(), o.first_called);
}

IGL_INLINE void build_lhs_numeric_multi(OPTICALData &o)
{
    result_file << "START BUILDING LEFT HAND SIDE USING NUMERIC TYPE MULTI-THREADED\n";
    igl::Timer t;
    t.start();
    compute_ix(o);
    compute_iy(o);
    compute_it(o);
    t.stop();
    write_to_file(result_file, "COMPUTE IX, IY, IT", t.getElapsedTimeInMicroSec(), o.first_called);
    if (o.first_called)
    {
        std::vector<Eigen::Triplet<ie::NumericType>> final_triplet;
        final_triplet.reserve(o.width * o.height * 3);
        o.datas.resize(2);
        o.datas[0].resize(o.width * o.height);
        o.datas[1].resize(o.width * o.height);
        int index = 0;
        ie::NumericType a2_numeric = ie::NumericType(o.alpha * o.alpha);
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, index, ie::NumericType(0, index) * ie::NumericType(0, index) + a2_numeric));
                index++;
            }
        }
        int index2 = index;
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, index, ie::NumericType(1, index - index2) * ie::NumericType(1, index - index2) + a2_numeric));
                index++;
            }
        }
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index2, index2 - o.height * o.width, ie::NumericType(0, index2 - o.height * o.width) * ie::NumericType(1, index2 - o.height * o.width)));
                index2++;
            }
        }
        Eigen::SparseMatrix<ie::NumericType> result_numeric;
        result_numeric.resize(o.width * o.height * 2, o.width * o.height * 2);
        result_numeric.setFromTriplets(final_triplet.begin(), final_triplet.end());
        result_numeric.makeCompressed();
        o.ex = ie::NumericExecutor(result_numeric, 0);
        o.result_vector.resize(result_numeric.nonZeros(), 0);
        // copy the outer index pointer and inner index pointer to o
        o.L_outer.resize(result_numeric.rows() + 1);
        o.L_inner.resize(result_numeric.nonZeros());
        o.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
        o.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
        ie::NumericType::clear_pool();
    }
    t.start();
    o.datas[0].assign(o.Ix.data(), o.Ix.data() + o.width * o.height);
    o.datas[1].assign(o.Iy.data(), o.Iy.data() + o.width * o.height);
    o.ex.ExecuteMulti(o.datas, o.result_vector);
    t.stop();
    write_to_file(result_file, "NUMERIC MULTI COMPUTE SPARSE MATRIX", t.getElapsedTimeInMicroSec(), o.first_called);
}

IGL_INLINE void build_lhs_numeric_multi_together(OPTICALData &o)
{
    result_file << "START BUILDING LEFT HAND SIDE USING NUMERIC TYPE BUILDING EVERYTHING TOGETHER\n";
    igl::Timer t;
    if (o.first_called)
    {
        typedef Eigen::Matrix<ie::NumericType, Eigen::Dynamic, Eigen::Dynamic> MatrixNum;
        MatrixNum image1_numeric, image2_numeric, ix_numeric, iy_numeric, it_numeric, u_bar_numeric, v_bar_numeric, u_numeric, v_numeric;
        image1_numeric.resize(o.height, o.width);
        image2_numeric.resize(o.height, o.width);
        ix_numeric.resize(o.height, o.width);
        iy_numeric.resize(o.height, o.width);
        it_numeric.resize(o.height, o.width);
        u_bar_numeric.resize(o.height, o.width);
        v_bar_numeric.resize(o.height, o.width);
        u_numeric.resize(o.height, o.width);
        v_numeric.resize(o.height, o.width);
        image1_numeric = ie::to_dense_numeric(o.image1, 0);
        image2_numeric = ie::to_dense_numeric(o.image2, 1);
        u_numeric = ie::to_dense_numeric(o.u, 2);
        v_numeric = ie::to_dense_numeric(o.v, 3);

        for (int w = 0; w < o.width - 1; w++)
        {
            for (int h = 0; h < o.height - 1; h++)
            {
                ix_numeric(h, w) = (image1_numeric(h, w + 1) - image1_numeric(h, w) + image1_numeric(h + 1, w + 1) - image1_numeric(h + 1, w) + image2_numeric(h, w + 1) - image2_numeric(h, w) + image2_numeric(h + 1, w + 1) - image2_numeric(h + 1, w)) / 4.0;
                iy_numeric(h, w) = (image1_numeric(h + 1, w) - image1_numeric(h, w) + image1_numeric(h + 1, w + 1) - image1_numeric(h, w + 1) + image2_numeric(h + 1, w) - image2_numeric(h, w) + image2_numeric(h + 1, w + 1) - image2_numeric(h, w + 1)) / 4.0;
            }
        }
        // set the boundary now
        for (int w = 0; w < o.width - 1; w++)
        {
            ix_numeric(o.height - 1, w) = 0.0;
            iy_numeric(o.height - 1, w) = 0.0;
        }
        for (int h = 0; h < o.height; h++)
        {
            ix_numeric(h, o.width - 1) = 0.0;
            iy_numeric(h, o.width - 1) = 0.0;
        }
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                it_numeric(h, w) = (return_matrix_value_numeric(image2_numeric, h, w) - return_matrix_value_numeric(image1_numeric, h, w) + return_matrix_value_numeric(image2_numeric, h + 1, w) - return_matrix_value_numeric(image1_numeric, h + 1, w) + return_matrix_value_numeric(image2_numeric, h, w + 1) - return_matrix_value_numeric(image1_numeric, h, w + 1) + return_matrix_value_numeric(image2_numeric, h + 1, w + 1) - return_matrix_value_numeric(image1_numeric, h + 1, w + 1)) / 4.0;
                u_bar_numeric(h, w) = (return_matrix_value_numeric(u_numeric, h + 1, w) + return_matrix_value_numeric(u_numeric, h, w + 1) + return_matrix_value_numeric(u_numeric, h - 1, w) + return_matrix_value_numeric(u_numeric, h, w - 1)) / 6.0 + (return_matrix_value_numeric(u_numeric, h - 1, w - 1) + return_matrix_value_numeric(u_numeric, h - 1, w + 1) + return_matrix_value_numeric(u_numeric, h + 1, w - 1) + return_matrix_value_numeric(u_numeric, h + 1, w + 1)) / 12.0;
                v_bar_numeric(h, w) = (return_matrix_value_numeric(v_numeric, h + 1, w) + return_matrix_value_numeric(v_numeric, h, w + 1) + return_matrix_value_numeric(v_numeric, h - 1, w) + return_matrix_value_numeric(v_numeric, h, w - 1)) / 6.0 + (return_matrix_value_numeric(v_numeric, h - 1, w - 1) + return_matrix_value_numeric(v_numeric, h - 1, w + 1) + return_matrix_value_numeric(v_numeric, h + 1, w - 1) + return_matrix_value_numeric(v_numeric, h + 1, w + 1)) / 12.0;
            }
        }

        std::vector<Eigen::Triplet<ie::NumericType>> final_triplet;
        std::vector<Eigen::Triplet<ie::NumericType>> rhs_triplet;
        final_triplet.reserve(o.width * o.height * 3);
        rhs_triplet.reserve(o.width * o.height * 2);
        o.datas.resize(4);
        o.datas[0].resize(o.width * o.height);
        o.datas[1].resize(o.width * o.height);
        o.datas[2].resize(o.width * o.height);
        o.datas[3].resize(o.width * o.height);
        int index = 0;
        ie::NumericType a2_numeric = ie::NumericType(o.alpha * o.alpha);
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, index, ix_numeric(h, w) * ix_numeric(h, w) + a2_numeric));
                index++;
            }
        }
        int index2 = index;
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, index, iy_numeric(h, w) * iy_numeric(h, w) + a2_numeric));
                index++;
            }
        }
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                final_triplet.push_back(Eigen::Triplet<ie::NumericType>(index2, index2 - o.height * o.width, ix_numeric(h, w) * iy_numeric(h, w)));
                index2++;
            }
        }
        index = 0;
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                rhs_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, 0, a2_numeric * u_bar_numeric(h, w) - ix_numeric(h, w) * it_numeric(h, w)));
                index++;
            }
        }
        for (int w = 0; w < o.width; w++)
        {
            for (int h = 0; h < o.height; h++)
            {
                rhs_triplet.push_back(Eigen::Triplet<ie::NumericType>(index, 0, a2_numeric * v_bar_numeric(h, w) - iy_numeric(h, w) * it_numeric(h, w)));
                index++;
            }
        }

        Eigen::SparseMatrix<ie::NumericType> result_numeric, rhs_numeric;
        result_numeric.resize(o.width * o.height * 2, o.width * o.height * 2);
        result_numeric.setFromTriplets(final_triplet.begin(), final_triplet.end());
        result_numeric.makeCompressed();
        rhs_numeric.resize(o.width * o.height * 2, 1);
        rhs_numeric.setFromTriplets(rhs_triplet.begin(), rhs_triplet.end());
        rhs_numeric.makeCompressed();
        o.ex = ie::NumericExecutor(result_numeric, 0);
        o.ex_rhs = ie::NumericExecutor(rhs_numeric, 0);
        o.result_vector.resize(result_numeric.nonZeros(), 0);
        o.rhs_vector.resize(rhs_numeric.nonZeros(), 0);
        // copy the outer index pointer and inner index pointer to o
        o.L_outer.resize(result_numeric.rows() + 1);
        o.L_inner.resize(result_numeric.nonZeros());
        o.L_outer.assign(result_numeric.outerIndexPtr(), result_numeric.outerIndexPtr() + result_numeric.rows() + 1);
        o.L_inner.assign(result_numeric.innerIndexPtr(), result_numeric.innerIndexPtr() + result_numeric.nonZeros());
        ie::NumericType::clear_pool();
    }
    t.start();
    o.datas[0].assign(o.image1.data(), o.image1.data() + o.width * o.height);
    o.datas[1].assign(o.image2.data(), o.image2.data() + o.width * o.height);
    o.datas[2].assign(o.u.data(), o.u.data() + o.width * o.height);
    o.datas[3].assign(o.v.data(), o.v.data() + o.width * o.height);
    o.ex.ExecuteMulti(o.datas, o.result_vector);
    o.ex_rhs.ExecuteMulti(o.datas, o.rhs_vector);
    t.stop();
    write_to_file(result_file, "NUMERIC MULTI COMPUTE EVERYTHING", t.getElapsedTimeInMicroSec(), o.first_called);
}

IGL_INLINE void build_rhs(OPTICALData &o)
{
    if (o.rhs.size() == 0)
    {
        o.rhs.resize(o.height * o.width * 2);
    }
    compute_ubar(o);
    compute_vbar(o);
    Eigen::MatrixXd rhs1 = std::pow(o.alpha, 2) * o.ubar - o.Ix.cwiseProduct(o.It);
    Eigen::MatrixXd rhs2 = std::pow(o.alpha, 2) * o.vbar - o.Iy.cwiseProduct(o.It);

    int index = 0;
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            o.rhs(index) = rhs1(h, w);
            index++;
        }
    }
    for (int w = 0; w < o.width; w++)
    {
        for (int h = 0; h < o.height; h++)
        {
            o.rhs(index) = rhs2(h, w);
            index++;
        }
    }
}

IGL_INLINE void solve_flow(OPTICALData &o)
{
    igl::Timer t;
    if (o.first_called)
    {
        result_file.open("result_opt.txt");
    }
    std::vector<double> solved(o.height * o.width * 2); // because we solve for u and v together
    switch (o.method)
    {
    case 0:
        build_lhs(o);
        t.start();
        build_rhs(o);
        t.stop();
        write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), o.first_called);
        break;
    case 3:
        build_lhs_numeric_multi(o);
        t.start();
        build_rhs(o);
        t.stop();
        write_to_file(result_file, "BUILDING RHS", t.getElapsedTimeInMicroSec(), o.first_called);
        break;
    case 4:
        build_lhs_numeric_multi_together(o);
        break;
    default:
        build_lhs(o);
    }

    if (o.first_called)
    {
        pardiso_init(o.pardiso_data);
        if (o.result_vector.size() == 0)
            pardiso_support_matrix(o.pardiso_data, o.lhs);
        else
            pardiso_support_matrix(o.pardiso_data, o.L_outer.data(), o.L_inner.data(), o.result_vector.data(), o.height * o.width * 2);
        pardiso_symbolic_factor(o.pardiso_data);
    }
    else
    {
        if (o.result_vector.size() == 0)
            pardiso_support_value(o.pardiso_data, o.lhs.valuePtr());
        else
            pardiso_support_value(o.pardiso_data, o.result_vector.data());
    }
    t.start();
    pardiso_numeric_factor(o.pardiso_data);
    t.stop();
    write_to_file(result_file, "MKL NUMERIC FACTOR", t.getElapsedTimeInMicroSec(), o.first_called);
    t.start();
    if (o.method <= 3)
    {
        pardiso_solve(o.pardiso_data, solved.data(), o.rhs.data());
    }
    else
    {
        pardiso_solve(o.pardiso_data, solved.data(), o.rhs_vector.data());
    }
    t.stop();
    write_to_file(result_file, "MKL SOLVE", t.getElapsedTimeInMicroSec(), o.first_called);

    t.start();
    o.u = Eigen::Map<Eigen::MatrixXd>(solved.data(), o.height, o.width);
    o.v = Eigen::Map<Eigen::MatrixXd>(solved.data() + o.height * o.width, o.height, o.width);
    t.stop();
    write_to_file(result_file, "PUTTING BACK TO U AND V", t.getElapsedTimeInMicroSec(), o.first_called);
    o.first_called = false;

    // compute_ix(o);
    // compute_iy(o);
    // compute_it(o);
    // // o.Ix *= -1;
    // // o.Iy *= -1;
    // // o.It *= -1;

    // Eigen::MatrixXd x2 = o.Ix.cwiseProduct(o.Ix);
    // Eigen::MatrixXd y2 = o.Iy.cwiseProduct(o.Iy);
    // for (int i = 0; i < 100; i++)
    // {
    //     std::cout << i << "\n";
    //     compute_ubar(o);
    //     compute_vbar(o);
    //     // std::cout<<o.Iy<<"\n";
    //     o.u = o.ubar - o.Ix.cwiseProduct((o.Ix.cwiseProduct(o.ubar) + o.Iy.cwiseProduct(o.vbar) + o.It).cwiseQuotient(x2 + y2 + Eigen::MatrixXd::Constant(o.height, o.width, 1.0)));
    //     o.v = o.vbar - o.Iy.cwiseProduct((o.Ix.cwiseProduct(o.ubar) + o.Iy.cwiseProduct(o.vbar) + o.It).cwiseQuotient(x2 + y2 + Eigen::MatrixXd::Constant(o.height, o.width, 1.0)));
    // }
}

} // namespace igl