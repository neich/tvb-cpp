//    Copyright 2020-2021 Ignacio Martín <ignacio.martin@udg.edu>
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//            http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#ifndef TVB_CPP_EQUATION_H
#define TVB_CPP_EQUATION_H

#include <string>
#include <unordered_map>

#include <tvb-cpp/exprtk.hpp>
#include <tvb-cpp/definitions.h>

namespace tvb {

    class Equation {
    public:
        typedef exprtk::symbol_table<double> SymbolTable;
        typedef exprtk::expression<double> Expression;
        typedef exprtk::parser<double> Parser;

    protected:
        SymbolTable m_symbol_table;
        Expression m_expression;

    public:
        Equation(const std::string &expression_string, const std::unordered_map<std::string, double> &constants) {
            for (auto &p: constants)
                m_symbol_table.add_constant(p.first, p.second);
            m_symbol_table.create_variable("X");
            m_expression.register_symbol_table(m_symbol_table);
            Parser parser;
            parser.compile(expression_string, m_expression);
        }

        virtual ~Equation() = default;

        double evaluate(const std::string& var, double value) const {
            m_symbol_table.get_variable(var)->ref() = value;
            return m_expression.value();
        }

        TArray1d evaluate(const std::string& var, const TArray1d& values) const {
            TArray1d result(values.size());
            for (unsigned i = 0; i < values.size(); ++i) {
                m_symbol_table.get_variable(var)->ref() = values[i];
                result[i] = m_expression.value();
            }
            return result;
        }


        double evaluate(const std::unordered_map<std::string, double> &parameters) {
            for (auto &p: parameters) {
                m_symbol_table.create_variable(p.first);
                m_symbol_table.get_variable(p.first)->ref() = p.second;
            }
            return m_expression.value();
        }

        double getVariableValue(const std::string& var) const {
            return m_symbol_table.get_variable(var)->ref();
        }
    };

    class HRFKernelEquation : public Equation {
    public:
        HRFKernelEquation(const std::string &expression_string,
                          const std::unordered_map<std::string, double> &parameters) :
                Equation(expression_string, parameters) {}
    };

    class FirstOrderVolterra : public HRFKernelEquation {
    public:
        FirstOrderVolterra() : HRFKernelEquation(
                "1/3.0 * exp(-0.5*(X / tau_s)) * (sin(sqrt(1.0/tau_f - 1.0/(4.0*tau_s*tau_s)) * X)) / (sqrt(1.0/tau_f - 1.0/(4.0*tau_s*tau_s)))",
                {{"tau_s", 0.8},
                 {"tau_f", 0.4},
                 {"k_1",   5.6},
                 {"V_0",   0.02}}) {}
    };
}


#endif //TVB_CPP_EQUATION_H
