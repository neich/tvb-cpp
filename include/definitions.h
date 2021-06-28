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

#define EIGEN_USE_MKL_ALL 1

#include <iostream>
#include <fstream>
#include <string>
#include <cstdarg>
#include <memory>
#include <unordered_map>

#include <Eigen/Dense>

namespace tvb {

    typedef typename std::complex<double> complexd;
    // Eigen::ColMajor is the default
    typedef typename Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixd;
    typedef typename Eigen::ArrayXXi Matrixi;
    typedef typename Eigen::Array<double, Eigen::Dynamic, 1, Eigen::ColMajor> Vectord;
    typedef typename Eigen::ArrayXi Vectori;
    typedef typename Eigen::MatrixXd AMatrixd;
    typedef typename Eigen::ArrayXd AVectord;
    typedef typename Eigen::Array<complexd, Eigen::Dynamic,1> AVectorc;
    typedef typename std::unordered_map<std::string, Matrixd> MatrixdMap;

//    typedef typename Eigen::MatrixXd Matrixd;
//    typedef typename Eigen::MatrixXi Matrixi;
//    typedef typename Eigen::VectorXd Vectord;
//    typedef typename Eigen::VectorXi Vectori;
//    typedef typename Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> AMatrixd;
//    typedef typename Eigen::ArrayXd AVectord;
//    typedef typename Eigen::Array<complexd, Eigen::Dynamic,1> AVectorc;
//    typedef typename std::unordered_map<std::string, Matrixd> MatrixdMap;


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
    Vectord operator/(const Vectord &v1, const Vectord &v2) {
        Vectord result(v1.size());
        for (unsigned i = 0; i < v1.size(); ++i)
            result(i) = v1(i) / v2(i);
        return result;
    }

    inline
    Vectord operator-(const double d, const Vectord &v) {
        Vectord result(v.size());
        for (unsigned i = 0; i < v.size(); ++i)
            result(i) = d - v(i);
        return result;
    }

    inline
    int index_circ(int index, int mod, int shift) {
        int result = (index % mod) + shift;
        return result >= 0 ? result : result + mod;
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
