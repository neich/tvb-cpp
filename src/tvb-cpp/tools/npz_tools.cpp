//
// Created by imartin on 07-Oct-20.
//

#include "npz_tools.h"

using namespace tvb;

template<typename Numeric>
TArray2d load_data_matrix(cnpy::NpyArray &w_npy) {
    unsigned rows = w_npy.shape[0];
    unsigned cols = w_npy.shape[1];

    auto* loaded_data = w_npy.data<Numeric>();

    TArray2d w(rows, cols);

    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < cols; ++c)
            w(r, c) = loaded_data[r + c * rows];

    return w;
}

template<typename Numeric>
TArray2d load_data_vector(cnpy::NpyArray &w_npy) {
    unsigned size = w_npy.shape[0];

    auto* loaded_data = w_npy.data<Numeric>();

    TArray1d w(size);

    for (unsigned i = 0; i < size; ++i)
            w[i] = loaded_data[i];

    return w;
}


TArray2d tvb::npz2Matrixd(const std::string& filename, const std::string& index) {
    cnpy::NpyArray w_npy = cnpy::npz_load(filename, index);
    assert(w_npy.shape.size() == 2);
    if (w_npy.word_size == sizeof(float)) return load_data_matrix<float>(w_npy);
    else if (w_npy.word_size == sizeof(double )) return load_data_matrix<double>(w_npy);
    else throw std::runtime_error(string_format("Unknown data format for file <%s>", filename.c_str()));
}

TArray2d tvb::npy2Matrixd(const std::string& filename) {
    cnpy::NpyArray w_npy = cnpy::npy_load(filename);
    assert(w_npy.shape.size() == 2);
    if (w_npy.word_size == sizeof(float)) return load_data_matrix<float>(w_npy);
    else if (w_npy.word_size == sizeof(double )) return load_data_matrix<double>(w_npy);
    else throw std::runtime_error(string_format("Unknown data format for file <%s>", filename.c_str()));
}

TArray1d tvb::npy2Vector(const std::string& filename) {
    cnpy::NpyArray w_npy = cnpy::npy_load(filename);
    assert(w_npy.shape.size() == 1);
    if (w_npy.word_size == sizeof(float)) return load_data_vector<float>(w_npy);
    else if (w_npy.word_size == sizeof(double )) return load_data_vector<double>(w_npy);
    else throw std::runtime_error(string_format("Unknown data format for file <%s>", filename.c_str()));
}



TArray2dMap tvb::npz2MatrixdMap(const std::string& filename) {
    cnpy::npz_t npy_map = cnpy::npz_load(filename);

    TArray2dMap result;
    for (auto& [key, w_npy]: npy_map) {
        assert(w_npy.shape.size() == 2);
        assert(w_npy.word_size == sizeof(Float));

        unsigned rows = w_npy.shape[0];
        unsigned cols = w_npy.shape[1];

        Float *loaded_data = w_npy.data<Float>();

        TArray2d w(rows, cols);

        for (unsigned r = 0; r < rows; ++r)
            for (unsigned c = 0; c < cols; ++c)
                w(r, c) = loaded_data[r + c * rows];

        result[key] = w;
    }

    return result;
}


std::vector<TArray2d> tvb::np2VecMatrixd(const std::string& filename, const std::string& index) {
    cnpy::NpyArray w_npy = index.size() == 0 ? cnpy::npy_load(filename) : cnpy::npz_load(filename, index);
    assert(w_npy.shape.size() == 3);
    assert(w_npy.word_size == sizeof(double));

    unsigned d0 = w_npy.shape[0];
    unsigned d1 = w_npy.shape[1];
    unsigned d2 = w_npy.shape[2];

    auto* loaded_data = w_npy.data<double>();

    std::vector<TArray2d> result(d0);
    std::fill(result.begin(), result.end(), tvb::TArray2d(d1, d2));
    for (unsigned i0 = 0; i0 < d0; ++i0)
        for (unsigned i1 = 0; i1 < d1; ++i1)
            for (unsigned i2 = 0; i2 < d2; ++i2)
                result[i0](i1, i2) = *(loaded_data++);

    return result;
}

std::vector<double> tvb::npz2VecDouble(const std::string& filename, const std::string& index) {
    cnpy::NpyArray w_npy = cnpy::npz_load(filename, index);
    assert(w_npy.shape.size() == 1);
    assert(w_npy.word_size == sizeof(double));

    double* loaded_data = w_npy.data<double>();

    std::vector<double> result;
    result.reserve(w_npy.shape[0]);
    for (unsigned i0 = 0; i0 < w_npy.shape[0]; ++i0) {
        result.push_back(loaded_data[i0]);
    }
    return result;
}

void tvb::TArray1d2npz(const TArray1d& vec, const std::string& filename, const std::string& index) {
    unsigned mat_size = vec.rows() * vec.cols();
    cnpy::npz_save<Float>(filename, index, vec.data(), {(size_t)vec.rows(), (size_t)vec.cols()});
}

void tvb::vecMatrixd2npz(const std::vector<TArray2d>& data, const std::string& filename, const std::string& index) {
    unsigned mat_size = data[0].rows()*data[0].cols();
    double* raw_data = new double[data.size()*mat_size];
    for (unsigned i = 0; i < data.size(); ++i)
        memcpy(&raw_data[i*mat_size], data[i].data(), mat_size*sizeof(double));
    cnpy::npz_save(filename, index, raw_data, {data.size(), (size_t)data[0].rows(), (size_t)data[0].cols()});
    delete[] raw_data;
}

void tvb::Matrixd2np(const TArray2d& data, const std::string& filename, const std::string& index) {
    unsigned mat_size = data.rows()*data.cols();
    // TODO: it assumes default ColMajor
    double* raw_data = new double[mat_size];
    for (int col = 0; col < data.cols(); ++col)
        for (int row = 0; row < data.rows(); ++row)
        raw_data[col*data.rows() + row] = data(row, col);
    if (index.size() > 0)
        cnpy::npz_save(filename, index, raw_data, { (size_t)data.rows(),  (size_t)data.cols()});
    else
        cnpy::npy_save(filename, raw_data, { (size_t)data.rows(),  (size_t)data.cols()});
    delete[] raw_data;
}


void tvb::vecDouble2npz(const std::vector<double>& data, const std::string& filename, const std::string& index) {
    cnpy::npz_save(filename, index, data.data(), {data.size()}, "a");
}

void tvb::MatrixdMap2npz(const std::string &filename, const TArray2dMap& mmap) {
    bool first = true;
    for (auto& [key, matrix]: mmap) {
        std::string mode = first ? "w" : "a";
        Eigen::Matrix<tvb::Float, -1, -1, Eigen::RowMajor> mrowmajor = matrix;
        cnpy::npz_save(filename, key, mrowmajor.data(), { (size_t)mrowmajor.rows(),  (size_t)mrowmajor.cols()}, mode);
        first = false;
    }
}


