//
// Created by imartin on 26-Oct-22.
//

#ifndef TVB_ROOT_CPP_PROJECTS_MODIFIER_H
#define TVB_ROOT_CPP_PROJECTS_MODIFIER_H


#include <tvb-cpp/simulator/model.h>

namespace tvb {

    class ParamModifier {
        int m_n;
    public:
        ParamModifier(int n): m_n(n) {}

        virtual TArray1d operator()(float t) = 0;
    };


    class ModelModifier {
    protected:
        Float m_t_start;
        Float m_t_end;
        Float m_dt;

        std::unordered_map<std::string, ParamModifier *> m_modifiers;

    public:
        ModelModifier() = default;

        virtual ~ModelModifier() = default;

        void configure(Float t_start, Float t_end, Float dt) {
            m_t_start = t_start;
            m_t_end = t_end;
            m_dt = dt;
        }

        void addModifier(const std::string& pname, ParamModifier* pm) {
            m_modifiers[pname] = pm;
        }

        void update(int step, Model *model) {
            float t = m_t_start + step * m_dt;
            for (auto [p, f]: m_modifiers) {
                model->set_param(p, (*f)(t));
            }
        };

    };

};
#endif //TVB_ROOT_CPP_PROJECTS_MODIFIER_H
