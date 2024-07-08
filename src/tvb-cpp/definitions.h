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
#include "tools/cnpy.h"

namespace tvb {

    typedef double Float;
    // typedef double Float;

    typedef typename std::complex<Float> complexd;
    // Eigen::ColMajor is the default
    typedef typename Eigen::Array<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> TArray2d;
    typedef typename std::unique_ptr<TArray2d> TArray2d_uptr;
    typedef typename std::shared_ptr<TArray2d> TArray2d_sptr;

    typedef typename Eigen::Array<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TArrayRM2d;
    typedef typename Eigen::Array<complexd, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> TArray2dc;
    typedef typename Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> TMatrix;
    typedef typename Eigen::ArrayXXi TArray2di;

    typedef typename Eigen::Array<Float, Eigen::Dynamic, 1, Eigen::ColMajor> TArray1d;
    typedef typename std::unique_ptr<TArray1d> TArray1d_uptr;
    typedef typename std::shared_ptr<TArray1d> TArray1d_sptr;

    typedef typename Eigen::Matrix<Float, Eigen::Dynamic, 1, Eigen::ColMajor> TVector;
    typedef typename Eigen::ArrayXi TArray1di;
    typedef typename Eigen::Array<complexd, Eigen::Dynamic, 1> TArray1dc;
    typedef typename std::map<std::string, TArray2d_sptr> TArray2dMap;
    typedef typename std::unordered_map<std::string, TArray2d> TArray2dUMap;

    typedef typename Eigen::Array<Float *, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrixdp;



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

//    inline
//    TArray1d operator/(const TArray1d &v1, const TArray1d &v2) {
//        TArray1d result(v1.size());
//        for (unsigned i = 0; i < v1.size(); ++i)
//            result(i) = v1(i) / v2(i);
//        return result;
//    }

    inline
    TArray1d operator-(const tvb::Float d, const TArray1d &v) {
        TArray1d result(v.size());
        for (unsigned i = 0; i < v.size(); ++i)
            result(i) = d - v(i);
        return result;
    }

    inline
    int index_circ(int index, int mod, int shift=0) {
        int result = (index % mod) + shift;
        return result >= 0 ? result : result + mod;
    }

    inline
    bool isnan(const TArray1d &vector) {
        for (auto const &v: vector)
            if (std::isnan(v))
                return true;
        return false;
    }

    inline
    bool isfinite(const TArray1d &vector) {
        for (auto const &v: vector)
            if (!std::isfinite(v))
                return false;
        return true;
    }

    inline
    bool isnan(const TArray2d &matrix) {
        for (auto const &v: matrix.reshaped())
            if (std::isnan(v))
                return true;
        return false;
    }

    inline
    bool isfinite(const TArray2d &matrix) {
        for (auto const &v: matrix.reshaped())
            if (!std::isfinite(v))
                return false;
        return true;
    }


    inline
    bool replace_nan(TArray1d &vector, Float value) {
        for (auto &v: vector)
            if (std::isnan(v))
                v = value;
        return false;
    }

    template<typename Numeric>
    std::vector<Numeric> range(Numeric start, Numeric finish, int intervals) {
        Numeric delta = (finish - start) / intervals;
        std::vector<Numeric> result(intervals + 1);
        Numeric current = start;
        for (int i = 0; i < intervals; ++i) {
            result[i] = current;
            current += delta;
        }
        result[intervals] = finish;
        return result;
    }

    template<typename Numeric>
    TArray1d nrange(Numeric start, Numeric delta, int n) {
        TArray1d result(n);
        Numeric current = start;
        for (int i = 0; i < n; ++i) {
            result[i] = current;
            current += delta;
        }
        return result;
    }

    template<typename Numeric>
    Eigen::Array<Numeric, Eigen::Dynamic,1> arange(Numeric start, Numeric end, Numeric step = 1) {
        assert(step != 0);
        Numeric d = abs(end - start);
        int N = int(d / step);
        if (step >= 1.0 && int(d) % int(step) != 0) N++;
        Eigen::Array<Numeric, Eigen::Dynamic,1> result(N);
        for (int i = 0; i < N; ++i, start+=step)
            result[i] = start;
        return result;
    }

    TArray2d getArray(cnpy::NpyArray &w_npy);
}


template<typename ... Args>
std::string string_format(const std::string &format, Args ... args) {
    size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
    if (size <= 0) { throw std::runtime_error("Error during formatting."); }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

template<typename ... Args>
void printff(const std::string &format, Args ... args) {
    if (sizeof...(Args) == 0)
        std::cout << format << std::flush;
    else
        std::cout << string_format(format, args...) << std::flush;
}


#define FE_0(WHAT)
#define FE_1(WHAT, X) WHAT(X)
#define FE_2(WHAT, X, ...) WHAT(X)FE_1(WHAT, __VA_ARGS__)
#define FE_3(WHAT, X, ...) WHAT(X)FE_2(WHAT, __VA_ARGS__)
#define FE_4(WHAT, X, ...) WHAT(X)FE_3(WHAT, __VA_ARGS__)
#define FE_5(WHAT, X, ...) WHAT(X)FE_4(WHAT, __VA_ARGS__)
#define FE_6(WHAT, X, ...) WHAT(X)FE_5(WHAT, __VA_ARGS__)
#define FE_7(WHAT, X, ...) WHAT(X)FE_6(WHAT, __VA_ARGS__)
#define FE_8(WHAT, X, ...) WHAT(X)FE_7(WHAT, __VA_ARGS__)
#define FE_9(WHAT, X, ...) WHAT(X)FE_8(WHAT, __VA_ARGS__)
#define FE_10(WHAT, X, ...) WHAT(X)FE_9(WHAT, __VA_ARGS__)
#define FE_11(WHAT, X, ...) WHAT(X)FE_10(WHAT, __VA_ARGS__)
#define FE_12(WHAT, X, ...) WHAT(X)FE_11(WHAT, __VA_ARGS__)
#define FE_13(WHAT, X, ...) WHAT(X)FE_12(WHAT, __VA_ARGS__)
#define FE_14(WHAT, X, ...) WHAT(X)FE_13(WHAT, __VA_ARGS__)
#define FE_15(WHAT, X, ...) WHAT(X)FE_14(WHAT, __VA_ARGS__)
#define FE_16(WHAT, X, ...) WHAT(X)FE_15(WHAT, __VA_ARGS__)
#define FE_17(WHAT, X, ...) WHAT(X)FE_16(WHAT, __VA_ARGS__)
#define FE_18(WHAT, X, ...) WHAT(X)FE_17(WHAT, __VA_ARGS__)
#define FE_19(WHAT, X, ...) WHAT(X)FE_18(WHAT, __VA_ARGS__)
#define FE_20(WHAT, X, ...) WHAT(X)FE_19(WHAT, __VA_ARGS__)
#define FE_21(WHAT, X, ...) WHAT(X)FE_20(WHAT, __VA_ARGS__)
#define FE_22(WHAT, X, ...) WHAT(X)FE_21(WHAT, __VA_ARGS__)
#define FE_23(WHAT, X, ...) WHAT(X)FE_22(WHAT, __VA_ARGS__)
#define FE_24(WHAT, X, ...) WHAT(X)FE_23(WHAT, __VA_ARGS__)
#define FE_25(WHAT, X, ...) WHAT(X)FE_24(WHAT, __VA_ARGS__)
#define FE_26(WHAT, X, ...) WHAT(X)FE_25(WHAT, __VA_ARGS__)
#define FE_27(WHAT, X, ...) WHAT(X)FE_26(WHAT, __VA_ARGS__)
#define FE_28(WHAT, X, ...) WHAT(X)FE_27(WHAT, __VA_ARGS__)
#define FE_29(WHAT, X, ...) WHAT(X)FE_28(WHAT, __VA_ARGS__)
#define FE_30(WHAT, X, ...) WHAT(X)FE_29(WHAT, __VA_ARGS__)
#define FE_31(WHAT, X, ...) WHAT(X)FE_30(WHAT, __VA_ARGS__)
#define FE_32(WHAT, X, ...) WHAT(X)FE_31(WHAT, __VA_ARGS__)
#define FE_33(WHAT, X, ...) WHAT(X)FE_32(WHAT, __VA_ARGS__)
#define FE_34(WHAT, X, ...) WHAT(X)FE_33(WHAT, __VA_ARGS__)

#define GET_MACRO(_0,_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,_12,_13,_14,_15,_16,_17,_18,_19,_20,_21,_22,_23,_24,_25,_26,_27,_28,_29,_30,_31,_32,_33,_34,NAME,...) NAME
#define FOR_EACH(action,...) \
      GET_MACRO(_0,__VA_ARGS__,FE_34,FE_33,FE_32,FE_31,FE_30,FE_29,FE_28,FE_27,FE_26,FE_25,FE_24,FE_23,FE_22,FE_21,FE_20,FE_19,FE_18,FE_17,FE_16,FE_15,FE_14,FE_13,FE_12,FE_11,FE_10,FE_9,FE_8,FE_7,FE_6,FE_5,FE_4,FE_3,FE_2,FE_1,FE_0)(action,__VA_ARGS__)

#define SETTER_FILL(field) if (#field == param) { this->field.fill(value); this->init_dependant(); return; }
#define SETTER_VALUE(field) if (#field == param) { if (value.size() == 1) { this->field.fill(value[0]); this->init_dependant(); } else { this->field=value; this->init_dependant(); } return; }
#define SETTER_SCALAR(field) if (#field == param) { this->field=value; this->init_dependant(); return; }
#define GETTER(field) if (#field == param) { return this->field; }
#define ADD_QUOTE(field) #field,

#define ADD_SETTER_FILL(...) FOR_EACH(SETTER_FILL, __VA_ARGS__)
#define ADD_SETTER_VALUE(...) FOR_EACH(SETTER_VALUE, __VA_ARGS__)
#define ADD_SETTER_SCALAR(...) FOR_EACH(SETTER_SCALAR, __VA_ARGS__)
#define ADD_GETTER(...) FOR_EACH(GETTER, __VA_ARGS__)

#define ADD_GETTERS_AND_SETTERS_SCALAR(...) \
    std::vector<std::string> get_param_list() const override {\
        return {FOR_EACH(ADD_QUOTE, __VA_ARGS__)};\
}                                           \
\
void set_param(const std::string& param, tvb::Float value) override {\
    FOR_EACH(SETTER_SCALAR, __VA_ARGS__)\
    throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));\
}                                           \
\
tvb::Float get_param(const std::string& param) const override {\
    FOR_EACH(GETTER, __VA_ARGS__)\
    throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));\
}



#endif //TVB_CPP_DEFINITIONS_H
