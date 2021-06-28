//
// Created by imartin on 07-Oct-20.
//

#include <tools/npz_tools.h>

using namespace tvb;

Matrixd tvb::npz2Matrixd(const std::string& filename, const std::string& index) {
    cnpy::NpyArray w_npy = cnpy::npz_load(filename, index);
    assert(w_npy.shape.size() == 2);
    assert(w_npy.word_size == sizeof(double));

    unsigned rows = w_npy.shape[0];
    unsigned cols = w_npy.shape[1];

    double* loaded_data = w_npy.data<double>();

    Matrixd w(rows, cols);

    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < cols; ++c)
            w(r, c) = loaded_data[r + c * rows];

    return w;
}

MatrixdMap tvb::npz2MatrixdMap(const std::string& filename) {
    cnpy::npz_t npy_map = cnpy::npz_load(filename);

    MatrixdMap result;
    for (auto& [key, w_npy]: npy_map) {
        assert(w_npy.shape.size() == 2);
        assert(w_npy.word_size == sizeof(double));

        unsigned rows = w_npy.shape[0];
        unsigned cols = w_npy.shape[1];

        double *loaded_data = w_npy.data<double>();

        Matrixd w(rows, cols);

        for (unsigned r = 0; r < rows; ++r)
            for (unsigned c = 0; c < cols; ++c)
                w(r, c) = loaded_data[r + c * rows];

        result[key] = w;
    }

    return result;
}


std::vector<Matrixd> tvb::npz2VecMatrixd(const std::string& filename, const std::string& index) {
    cnpy::NpyArray w_npy = cnpy::npz_load(filename, index);
    assert(w_npy.shape.size() == 3);
    assert(w_npy.word_size == sizeof(double));

    unsigned d0 = w_npy.shape[0];
    unsigned d1 = w_npy.shape[1];
    unsigned d2 = w_npy.shape[2];

    auto* loaded_data = w_npy.data<double>();

    std::vector<Matrixd> result(d0);
    std::fill(result.begin(), result.end(), tvb::Matrixd(d1, d2));
    for (unsigned i2 = 0; i2 < d2; ++i2)
        for (unsigned i1 = 0; i1 < d1; ++i1)
                for (unsigned i0 = 0; i0 < d0; ++i0)
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


void tvb::vecMatrixd2npz(const std::vector<Matrixd>& data, const std::string& filename, const std::string& index) {
    unsigned mat_size = data[0].rows()*data[0].cols();
    double* raw_data = new double[data.size()*mat_size];
    for (unsigned i = 0; i < data.size(); ++i)
        memcpy(&raw_data[i*mat_size], data[i].data(), mat_size*sizeof(double));
    cnpy::npz_save(filename, index, raw_data, {data.size(), (size_t)data[0].rows(), (size_t)data[0].cols()}, "a");
    delete[] raw_data;
}

void tvb::Matrixd2npz(const Matrixd& data, const std::string& filename, const std::string& index) {
    unsigned mat_size = data.rows()*data.cols();
    // TODO: it assumes default ColMajor
    double* raw_data = new double[mat_size];
    for (int col = 0; col < data.cols(); ++col)
        for (int row = 0; row < data.rows(); ++row)
        raw_data[col*data.rows() + row] = data(row, col);
    cnpy::npz_save(filename, index, raw_data, { (size_t)data.rows(),  (size_t)data.cols()}, "a");
    delete[] raw_data;
}


void tvb::vecDouble2npz(const std::vector<double>& data, const std::string& filename, const std::string& index) {
    cnpy::npz_save(filename, index, data.data(), {data.size()}, "a");
}

void tvb::MatrixdMap2npz(const std::string &filename, const MatrixdMap& mmap) {
    bool first = true;
    for (auto& [key, matrix]: mmap) {
        std::string mode = first ? "w" : "a";
        Eigen::Matrix<double, -1, -1, Eigen::RowMajor> mrowmajor = matrix;
        cnpy::npz_save(filename, key, mrowmajor.data(), { (size_t)mrowmajor.rows(),  (size_t)mrowmajor.cols()}, mode);
        first = false;
    }
}


