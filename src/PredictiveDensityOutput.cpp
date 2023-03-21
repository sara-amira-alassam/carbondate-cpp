//
// Created by Sara Al-Assam on 20/03/2023.
//
#include "PredictiveDensityOutput.h"

std::vector<std::string> PredictiveDensityOutput::get_output_lines() {
    std::vector<std::string> output_lines = DensityOutput::get_output_lines();
    output_lines.push_back(_output_var + ".name=\"NP Model;\"\n");
    output_lines.push_back(_output_var + ".op=\"NP_Model\";\n");
    output_lines.push_back(_output_prefix + ".probNorm*=\"" + std::to_string(_n_obs) + "\";\n");
    return output_lines;
}

PredictiveDensityOutput::PredictiveDensityOutput(int n_obs, int offset, double resolution)
    : _n_obs(n_obs), DensityOutput(offset + 1, resolution) {
    _output_var = "ocd[" + std::to_string(_index + 1) + "]";
    _output_prefix = _output_var + ".posterior";
}

void PredictiveDensityOutput::set_confidence_intervals(
        const std::vector<double>& ci_lower,
        const std::vector<double>& ci_upper) {
    _ci_lower.resize(ci_lower.size());
    _ci_upper.resize(ci_upper.size());

    for (int i = 0; i < ci_lower.size(); i++) {
        _ci_lower[i] = ci_lower[i]/_prob_max;
        _ci_upper[i] = ci_upper[i]/_prob_max;
    }
}
