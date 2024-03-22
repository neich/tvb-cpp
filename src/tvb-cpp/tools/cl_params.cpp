//
// Created by natxm on 21/03/2024.
//

#include "cl_params.h"

#include <cassert>

#include <tvb-cpp/simulator/factory.h>
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/tools/npz_tools.h>

struct SweepParam {
    std::string name;
    std::vector<tvb::Float> values;
};

void tvb::CLParser::load_data() {
    tvb::TArray2d C;

    std::string file(vm[o_sc_matrix].as<std::string>());
    if (file.ends_with(".csv"))
        C = tvb::csv_load(file);
    else if (file.ends_with(".npz"))
        C = tvb::npz2Matrixd(file, "SC");
    else if (file.ends_with(".npy"))
        C = tvb::npy2Matrixd(file);
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file.c_str()));

    assert(("Connectivity matrix has to be square!", C.rows() == C.cols()));

    sc_matrix.reset(new tvb::TArray2d(C));

    n_rois = sc_matrix->rows();

    if (vm.contains(o_lenghts_matrix)) {
        std::string file_weights(vm[o_lenghts_matrix].as<std::string>());
        if (file_weights.ends_with(".csv"))
            C = tvb::csv_load(file_weights);
        else if (file_weights.ends_with(".npz"))
            C = tvb::npz2Matrixd(file_weights, "SC");
        else if (file_weights.ends_with(".npy"))
            C = tvb::npy2Matrixd(file_weights);
        else
            throw std::runtime_error(string_format("Unknown file extension for: %s", file_weights.c_str()));

        assert(("Lengths matrix has to be the same size than connectivity matrix", C.rows() == sc_matrix->rows() && C.cols() == sc_matrix->cols()));

        length_matrix.reset(new tvb::TArray2d(C));
    }
    else {
        length_matrix.reset(new TArray2d(n_rois, n_rois));
        // Make lengths very small so speed is equivalent to infinity and there is no delays
        length_matrix->setConstant(1e-10);
    }
}


std::vector<tvb::ParameterSet> tvb::CLParser::get_parameter_combinations() const {
    std::vector<SweepParam> params;
    if (vm.count(o_param)) {
        for (auto &s: vm[o_param].as<std::vector<std::string>>()) {
            if (std::isalpha(s[0])) {
                if (!params.empty() && !(params.back().values.size() == 1 || params.back().values.size() == 3))
                    throw std::runtime_error(string_format("Malformed parameter %s\n", s.c_str()));
                params.emplace_back();
                params.back().name = s;
            } else {
                try {
                    tvb::Float value = std::stof(s);
                    params.back().values.push_back(value);
                } catch (...) {
                    throw std::runtime_error(string_format("Syntax error in parameter value: <%s>", s.c_str()));
                }
            }
        }
        if (params.back().values.size() != 1 && params.back().values.size() != 3)
            throw std::runtime_error(string_format("Malformed parameter <%s>\n", params.back().name.c_str()));
    }

    std::vector<ParameterSet> param_combs(1);
    for (auto const &p: params) {
        if (p.values.size() == 1) {
            for (auto &pc: param_combs)
                pc.emplace_back(p.name, p.values[0]);
        } else {
            std::vector<ParameterSet> new_param_combs;
            for (auto v: tvb::range(p.values[0], p.values[1], p.values[2])) {
                for (auto const &pc: param_combs) {
                    new_param_combs.push_back(pc);
                    new_param_combs.back().emplace_back(p.name, v);
                }
            }
            param_combs = new_param_combs;
        }
    }

    return param_combs;
}


bpo::variables_map tvb::CLParser::init(bpo::options_description &desc, int argc, const char **argv) {
    desc.add_options()
            ("help,h", "Help screen")
            (o_sc_matrix.c_str(), bpo::value<std::string>()->required(), "Structural connectivity matrix")
            (o_param.c_str(), bpo::value<std::vector<std::string>>()->multitoken(), "Parameters to sweep")
            (o_model.c_str(), bpo::value<std::string>()->required(), "Whole brain model to use")
            (o_noise.c_str(), bpo::value<std::vector<tvb::Float>>()->multitoken(), "Vector with noise sigmas for each state variable");

    store(bpo::command_line_parser(argc, argv)
                  .options(desc)
                  .style(bpo::command_line_style::unix_style ^ bpo::command_line_style::allow_short)
                  .run(), vm);

    if (vm.count("help")) {
        std::cout << desc << '\n';
        exit(0);
    }

    notify(vm);

    load_data();

    model.reset(tvb::Factory::new_model(vm[o_model].as<std::string>(), n_rois));

    if (vm.contains(o_noise)) {
        auto v = vm[o_noise].as<std::vector<tvb::Float>>();
        noise = Eigen::Map<TArray1d, Eigen::Unaligned>(v.data(), v.size());
    }

    return vm;
}

tvb::Float tvb::CLParser::init_from_parameters(tvb::Model *model, const tvb::ParameterSet &pset) const {
    tvb::Float G = 1.0;
    for (auto const &p: pset)
        if (std::isalpha(p.name[0])) {
            if (p.name != scale_param)
                model->set_param(p.name, p.value);
            else
                G = p.value;
        }
    model->init_dependant();
    return G;
}

