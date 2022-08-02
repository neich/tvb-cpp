//    Copyright 2020-2021 Ignacio Martín <ignacio.martin@udg.edu>
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//            http://www.apache.org/licenses/LICENSE-2.0
//SI fue
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#ifndef TVB_CPP_DEFINITIONS_H
#define TVB_CPP_DEFINITIONS_H


#include <iostream>
#include <fstream>
#include <string>
#include <cstdarg>
#include <memory>
#include <unordered_map>

#include <Eigen/Dense>

namespace tvb {

    typedef float Float;
    // typedef double Float;

    typedef typename std::complex<Float> complexd;
    // Eigen::ColMajor is the default
    typedef typename Eigen::Array<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> TArray2d;
    typedef typename Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> TMatrix;
    typedef typename Eigen::ArrayXXi TArray2di;
    typedef typename Eigen::Array<Float, Eigen::Dynamic, 1, Eigen::ColMajor> TArray1d;
    typedef typename Eigen::Matrix<Float, Eigen::Dynamic, 1, Eigen::ColMajor> TVector;
    typedef typename Eigen::ArrayXi TArray1di;
    typedef typename Eigen::Array<complexd, Eigen::Dynamic,1> TArray1dc;
    typedef typename std::unordered_map<std::string, TArray2d> TArray2dMap;

    typedef typename Eigen::Array<Float*, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixdp;

//    typedef typename Eigen::MatrixXd TArray2d;
//    typedef typename Eigen::MatrixXi TArray2di;
//    typedef typename Eigen::VectorXd TArray1d;
//    typedef typename Eigen::VectorXi TArray1di;
//    typedef typename Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> AMatrixd;
//    typedef typename Eigen::ArrayXd AVectord;
//    typedef typename Eigen::Array<complexd, Eigen::Dynamic,1> TArray1dc;
//    typedef typename std::unordered_map<std::string, TArray2d> TArray2dMap;


    template<class M, class G>
    void generate(M &m, G g) {
        for (unsigned i = 0; i < m.rows(); ++i)
            for (unsigned j = 0; j < m.cols(); ++j)
                m(i, j) = g();
    }

    template<class M1, class M2, class G>
    void transform(const M1 &m_from, M2 &m_to, G g) {
        for (unsigned i = 0; i < m_from.rows(); ++i)
            for (unsigned j = 0; j < m_from.cols(); ++j) {
                int v = g(m_from(i, j));
                m_to(i, j) = v;
            }
    }

    template<class M>
    void fill(M &m, typename M::Scalar v) {
        for (unsigned i = 0; i < m.rows(); ++i)
            for (unsigned j = 0; j < m.cols(); ++j)
                m(i, j) = v;
    }

    template<class M>
    typename M::Scalar max(const M &m) {
        typename M::Scalar max = m(0, 0);
        for (unsigned i = 0; i < m.rows(); ++i)
            for (unsigned j = 0; j < m.cols(); ++j)
                if (m(i, j) > max)
                    max = m(i, j);
        return max;
    }

    inline
    TArray1d operator/(const TArray1d &v1, const TArray1d &v2) {
        TArray1d result(v1.size());
        for (unsigned i = 0; i < v1.size(); ++i)
            result(i) = v1(i) / v2(i);
        return result;
    }

    inline
    TArray1d operator-(const double d, const TArray1d &v) {
        TArray1d result(v.size());
        for (unsigned i = 0; i < v.size(); ++i)
            result(i) = d - v(i);
        return result;
    }

    inline
    int index_circ(int index, int mod, int shift) {
        int result = (index % mod) + shift;
        return result >= 0 ? result : result + mod;
    }

    inline
    bool isnan(const TArray1d& vector) {
        for (auto v: vector)
            if (std::isnan(v))
                return true;
        return false;
    }

    inline
    bool isnan(const TArray2d& matrix) {
        for (auto v: matrix.reshaped())
            if (std::isnan(v))
                return true;
        return false;
    }

    inline
    bool replace_nan(TArray1d& vector, Float value) {
        for (auto &v: vector)
            if (std::isnan(v))
                v = value;
        return false;
    }

    template<typename Numeric>
    std::vector<Numeric> range(Numeric start, Numeric finish, int intervals) {
        Numeric delta = (finish - start) / intervals;
        std::vector<Numeric> result(intervals+1);
        Numeric current = start;
        for (int i = 0; i < intervals; ++i) {
            result[i] = current;
            current += delta;
        }
        result[intervals] = finish;
        return result;
    }


}


template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf( new char[ size ] );
    snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

#endif //TVB_CPP_DEFINITIONS_H
