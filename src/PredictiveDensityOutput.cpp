#include <iostream>
#include <utility>
#include "helpers.h"
#include "PredictiveDensityOutput.h"

PredictiveDensityOutput::PredictiveDensityOutput(
        int n_obs, int offset, double resolution, std::string name,
        const std::vector<double>& cal_age_AD,
        const std::vector<double>& mean_density,
        const std::vector<double>& ci_lower,
        const std::vector<double>& ci_upper)
        : _n_obs(n_obs), _name(std::move(name)), DensityOutput(offset + 1, resolution) {

    if (cal_age_AD[1] - cal_age_AD[0] != resolution) {
        // We don't expect this to happen, but best to double-check
        std::cerr << "Resolution should be " + std::to_string(resolution);
        std::cerr << " but is " + std::to_string(cal_age_AD[1] - cal_age_AD[0]);
        exit(1);
    }

    set_probability(mean_density);
    set_confidence_intervals(ci_lower, ci_upper);

    _start_calAD = cal_age_AD[0];
    _mean_calAD = mean(cal_age_AD, mean_density);
    _median_calAD = median(cal_age_AD, mean_density);
    _sigma_calAD = sigma(cal_age_AD, mean_density, _mean_calAD);
}

std::vector<std::string> PredictiveDensityOutput::get_output_lines() {
    std::string model_line;
    double upper_lim = _start_calAD + _resolution * (double) _probability.size();
    model_line = "model.element[" + std::to_string(_index) + "]= ";
    model_line += R"({op: "NP_Model", type: "date", name: ")" + _name + "\", ";
    model_line += "pos:" + std::to_string(_index) + ", timepos: " + std::to_string(_index) + ", ";
    model_line += "lower: " + to_string(_start_calAD) + ", upper:" + to_string(upper_lim) + "};\n";

    std::vector<std::string> output_lines = DensityOutput::get_output_lines();
    output_lines.emplace_back(model_line);
    output_lines.push_back(_output_var + ".name=\"" + _name + "\"\n");
    output_lines.push_back(_output_var + ".op=\"NP_Model\";\n");
    output_lines.push_back(_output_prefix + ".probNorm*=\"" + std::to_string(_n_obs) + "\";\n");
    return output_lines;
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
