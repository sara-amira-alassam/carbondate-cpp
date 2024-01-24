/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include "PredictiveDensityOutput.h"

/*
 * Creates an object suitable for printing out the predictive calendar age density for all
 * determinations, taking the following arguments:
 * - n_obs:  The number of observations (this is used for normalization)
 * - offset: The offset to use for indexing the total model output
 *           (e.g. if other models have been specified in the input file in addition to this one)
 * - resolution: The resolution to use for determining the probability curve
 * - model_name: The name given to the model by the user
 * - cal_age_AD: A vector of the calendar ages to use
 * - mean_density: The mean of the sampled densities
 * - ci_lower, ci_upper: The 1-sigma confidence intervals of the sampled densities
 */
PredictiveDensityOutput::PredictiveDensityOutput(
        int n_obs, int offset, double resolution, std::string model_name,
        const std::vector<double>& cal_age_AD,
        const std::vector<double>& mean_density,
        const std::vector<double>& ci_lower,
        const std::vector<double>& ci_upper)
        : _n_obs(n_obs), _name(std::move(model_name)), DensityOutput(offset, resolution) {

    if (cal_age_AD[1] - cal_age_AD[0] != resolution) {
        // We don't expect this to happen, but best to double-check
        throw DensityOutputException(
                "Resolution should be " + std::to_string(resolution) +
                " but is " + std::to_string(cal_age_AD[1] - cal_age_AD[0]));
    }

    _set_probability(mean_density);
    set_confidence_intervals(ci_lower, ci_upper);

    _start_calAD = cal_age_AD[0];

    std::vector<double> normalised_density(mean_density.begin(), mean_density.end());
    double prob_total = 0.;
    for (double i : normalised_density) prob_total += i;
    for (double & i : normalised_density) i /= prob_total;
    _mean_calAD = mean(cal_age_AD, normalised_density);
    _median_calAD = median(cal_age_AD, normalised_density);
    _sigma_calAD = sigma(cal_age_AD, normalised_density, _mean_calAD);
}

void PredictiveDensityOutput::_generate_output_lines() {
    DensityOutput::_generate_output_lines();
    std::string model_line;
    double upper_lim = _start_calAD + _resolution * (double) _probability.size();
    model_line = "model.element[" + std::to_string(_index) + "]= ";
    model_line += R"({op: "NP_Model", type: "date", name: ")" + _name + "\", ";
    model_line += "pos:" + std::to_string(_index) + ", timepos: " + std::to_string(_index) + ", ";
    model_line += "lower: " + to_string(_start_calAD, 6) + ", upper:" + to_string(upper_lim, 6) + "};\n";

    std::string ci_line = "model.element[" + std::to_string(_index) + "].kde_mean={";
    ci_line += "probNorm:" + to_string(_prob_norm * _n_obs, 6) + ", ";
    ci_line += "start:" + to_string(_start_calAD, 6) + ", ";
    ci_line += "resolution:" + to_string(_resolution, 6) + ", ";
    ci_line += "prob_sigma:[";
    for (int i = 0; i < _ci_lower.size(); i++) {
        ci_line += "[" + to_string((_ci_upper[i] + _ci_lower[i]) / 2., 6);
        ci_line += "," + to_string((_ci_upper[i] - _ci_lower[i]) / 2., 6) + "]";
        if (i != _ci_lower.size() - 1) ci_line += ",";
    }
    ci_line += "]};\n";

    _js_output_lines.emplace_back(model_line);
    _js_output_lines.emplace_back(ci_line);
    _js_output_lines.push_back(_variable_line("name", _name));
    _js_output_lines.push_back(_variable_line("op", "NP_Model"));
    _js_output_lines.push_back(_output_line("probNorm", _prob_norm * _n_obs));
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

std::string PredictiveDensityOutput::_range_lines(int &comment_index) {
    return _output_prefix + ".range=[];\n";
}
