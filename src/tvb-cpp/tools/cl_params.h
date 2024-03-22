//
// Created by natxm on 21/03/2024.
//

#include <tvb-cpp/simulator/model.h>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

namespace tvb {

    struct Parameter {
        std::string name;
        tvb::Float value;

        Parameter(std::string name, tvb::Float value) : name(std::move(name)), value(value) {}

        bool operator==(const std::string &v) const { return name == v; }
    };

    typedef std::vector<Parameter> ParameterSet;

    class CLParser {
        std::string o_model = "model";
        std::string o_param = "param";
        std::string o_sc_matrix = "sc-matrix";
        std::string o_lenghts_matrix = "lengths-matrix";
        std::string o_noise = "noise";
        std::string scale_param = "G";
        tvb::Float G = 1.0;

        bpo::variables_map vm;

        int n_rois = 0;

        std::shared_ptr<tvb::TArray2d> sc_matrix;
        std::shared_ptr<tvb::TArray2d> length_matrix;
        std::shared_ptr<tvb::Model> model;
        tvb::TArray1d noise;

        void load_data();

    public:

        bpo::variables_map init(bpo::options_description& desc, int argc, const char **argv);

        bpo::variables_map& get_variables_map() {
            return vm;
        }

        const tvb::TArray1d& get_noise() const {
            return noise;
        }

        int get_n_rois() const {
            return n_rois;
        }

        std::shared_ptr<tvb::TArray2d> get_sc_matrix() const {
            return sc_matrix;
        }

        std::shared_ptr<tvb::TArray2d> get_length_matrix() const {
            return length_matrix;
        }


        [[nodiscard]] std::vector<ParameterSet> get_parameter_combinations() const;

        tvb::Float init_from_parameters(tvb::Model *model, const tvb::ParameterSet &pset) const;

        template<typename T>
        T get_option(std::string option) {
            return vm[option].as<T>();
        }
    };


}