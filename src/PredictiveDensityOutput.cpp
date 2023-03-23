#include <iostream>
#include <utility>
#include "helpers.h"
#include "PredictiveDensityOutput.h"

// Creates an object suitable for printing out the predictive calendar age density for all
// determinations, taking the following arguments:
// * n_obs:  The number of observations (this is used for normalization)
// * offset: The offset to use for indexing the total model output
//           (e.g. if other models have been specified in the input file in addition to this one)
// * resolution: The resolution to use for detmining the probability curve
// * mean_density: The mean of the sampled densities
// * ci_lower, ci_upper: The 1-sigma confidence intervals of the sampled densities
PredictiveDensityOutput::PredictiveDensityOutput(
        int n_obs, int offset, double resolution, std::string name,
        const std::vector<double>& cal_age_AD,
        const std::vector<double>& mean_density,
        const std::vector<double>& ci_lower,
        const std::vector<double>& ci_upper)
        : _n_obs(n_obs), _name(std::move(name)), DensityOutput(offset, resolution) {

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
    std::vector<std::string> output_lines = DensityOutput::get_output_lines();
    std::string model_line;
    double upper_lim = _start_calAD + _resolution * (double) _probability.size();
    model_line = "model.element[" + std::to_string(_index) + "]= ";
    model_line += R"({op: "NP_Model", type: "date", name: ")" + _name + "\", ";
    model_line += "pos:" + std::to_string(_index) + ", timepos: " + std::to_string(_index) + ", ";
    model_line += "lower: " + to_string(_start_calAD) + ", upper:" + to_string(upper_lim) + "};\n";

    std::string ci_line = "model.element[" + std::to_string(_index) + "].kde_mean={";
    ci_line += "probNorm:" + to_string(_prob_norm * _n_obs) + ", ";
    ci_line += "start:" + to_string(_start_calAD) + ", ";
    ci_line += "resolution:" + to_string(_resolution) + ", ";
    ci_line += "prob_sigma:[";
    for (int i = 0; i < _ci_lower.size(); i++) {
        ci_line += "[" + to_string(_ci_upper[i]) + "," + to_string(_ci_lower[i]) + "]";
        if (i != _ci_lower.size() - 1) ci_line += ",";
    }
    ci_line += "]};\n";

    output_lines.emplace_back(model_line);
    output_lines.emplace_back(ci_line);
    output_lines.push_back(variable_line("name", _name));
    output_lines.push_back(variable_line("op", "NP_Model"));
    output_lines.push_back(output_line("probNorm", _prob_norm * _n_obs));
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
