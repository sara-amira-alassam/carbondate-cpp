//
// Created by Sara Al-Assam on 21/02/2023.
//
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include "helpers.h"
#include "DensityOutput.h"

DensityOutput::DensityOutput(int index, double resolution)
    : _index(index), _resolution(resolution) {
    _output_var = "ocd[" + std::to_string(_index) + "]";
    _output_prefix = _output_var + ".posterior";
}

void DensityOutput::print(const std::string& file_prefix) {
    std::ofstream output_file;
    output_file.open("../output/" + file_prefix + ".js", std::ios_base::app);

    for (const std::string& output_line : get_output_lines()) {
        output_file << output_line;
    }
    output_file.close();
}

std::vector<std::string> DensityOutput::get_output_lines() {
    std::vector<std::string> output_lines;
    int comment_index = 0;

    output_lines.push_back("if(!" + _output_var + "){" + _output_var + "={};}\n");
    output_lines.push_back(variable_line("ref", "TODO: Add custom ref"));
    output_lines.push_back(_output_prefix + "={};\n");
    output_lines.push_back(comment_line("Posterior ", comment_index));
    output_lines.push_back(range_lines(comment_index));
    output_lines.push_back(output_line("mean", _mean_calAD));
    output_lines.push_back(output_line("sigma", _sigma_calAD));
    output_lines.push_back(output_line("median", _median_calAD));
    output_lines.push_back(output_line("probNorm", _prob_norm));
    output_lines.push_back(output_line("start", _start_calAD));
    output_lines.push_back(output_line("resolution", _resolution));
    output_lines.push_back(output_line("prob", _probability));

    return output_lines;
}

std::string DensityOutput::comment_line(const std::string& comment, int& comment_index) {
    std::string line;
    if (comment_index == 0) {
        line = _output_prefix + ".comment=[];\n";
    }
    line += _output_prefix + ".comment[" + std::to_string(comment_index++) + "]=\"" +  comment + "\";\n";
    return line;
}

std::string DensityOutput::variable_line(const std::string& var_name, const std::string& var) {
    return _output_var + "." + var_name + "=\"" +  var + "\";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, double var) {
    return _output_prefix + "." + var_name + "=" + to_string(var, 6) + ";\n";
}

std::string DensityOutput::output_line(
        const std::string& var_name, const std::vector<double>& var) {
    std::string output_line = _output_prefix + "." + var_name + "=[";
    for (int i = 0; i < var.size() - 1; i++) output_line += std::to_string(var[i]) + ", ";
    output_line += std::to_string(var[var.size() - 1]) + "];\n";
    return output_line;
}

void DensityOutput::set_probability(const std::vector<double>& probability) {
    double prob_total = 0.;
    _probability.resize(probability.size());
    for (int i = 0; i < probability.size(); i++) {
        _probability[i] = probability[i];
        prob_total += probability[i];
        if (probability[i] > _prob_max) _prob_max = probability[i];
    }
    for (int i = 0; i < probability.size(); i++) _probability[i] /= _prob_max;
    _prob_norm = _prob_max / (prob_total * _resolution);
}

// Returns the area under the probability curve if we ignore all values below the cut-off
// Also populates the vector ranges, where each entry contains
// [start_calAD, end_calAD, probability within this range]
double DensityOutput::find_probability_and_ranges_for_cut_off(
        double cut_off, std::vector<std::vector<double>> &ranges) {
    ranges.clear();
    double y1, y2, dx, x_intercept_1, x_intercept_2, res = _resolution;
    double range_probability = 0, total_probability = 0;
    const double min_prob = 0.005; // Don't bother to store ranges with probability less than this
    for (int i = 0; i < _probability.size() - 1; i++) {
        y1 = _probability[i];
        y2 = _probability[i + 1];
        if (y1 <= cut_off and y2 > cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_1 = _start_calAD + i * res + dx;
            range_probability = (cut_off + y2) * (res - dx) / 2.;
        } else if (y1 > cut_off and y2 <= cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_2 = _start_calAD + i * res + dx;
            range_probability += (y1 + cut_off) * dx / 2.;
            range_probability *= _prob_norm;
            if (range_probability > min_prob) {
                ranges.push_back(
                        std::vector<double> {x_intercept_1, x_intercept_2, range_probability});
                total_probability += range_probability;
            }
            range_probability = 0;
        } else if (y1 > cut_off and y2 > cut_off) {
            range_probability += (y1 + y2) * res / 2.;
        }
    }
    for (std::vector<double> & range : ranges) range[2] /= total_probability;
    return total_probability;
}

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges_by_intercepts(double probability) {
    std::vector<std::vector<double>> ranges;

    // Use bisection method to find the closest probability cut-off to give the desired probability
    // a and b are the upper and lower points of the section - we know the smoothed probability has a max of 1
    double a = 0., b = 1.;
    double p; // p is the midpoint between a and b
    double current_probability;
    const int max_iter = 1000;
    for (int i = 0; i < max_iter; i++) {
        p = (a + b) / 2.;
        current_probability = find_probability_and_ranges_for_cut_off(p, ranges);
        if (abs(current_probability - probability) < 1e-4) {
            break;
        }
        if (current_probability < probability) {
            b = p;
        } else {
            a = p;
        }
    }
    return ranges;
}

std::string DensityOutput::range_lines(int &comment_index) {
    return {};
}

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges(double probability) {
    double resolution = _resolution, sum_prob = 1 / (_prob_norm * _resolution);
    std::vector<double> probability_interp(_probability.begin(), _probability.end());
    int n = (int) _probability.size();

    // TODO: Really we want to create a smoothed density from the posterior values but here we
    // do a much simpler approximation of the histogram values
    if (_resolution > 1.0) {
        int scale_factor = floor(_resolution / 0.5);
        resolution = _resolution / scale_factor;
        sum_prob = 0.;
        n = ((int) _probability.size() - 1) * scale_factor + 1;
        probability_interp.resize(n, 0);
        double cal_age_interp, cal_age_beg, cal_age_end;
        for (int i = 0; i < _probability.size() - 1; i++) {
            cal_age_beg = _start_calAD + i * _resolution;
            cal_age_end = cal_age_beg + _resolution;
            for (int j = 0; j < scale_factor; j++) {
                cal_age_interp = cal_age_beg + j * resolution;
                probability_interp[i * scale_factor + j] = interpolate_linear(
                    cal_age_interp, cal_age_beg, cal_age_end, _probability[i], _probability[i + 1]);
                sum_prob += probability_interp[i * scale_factor + j];
            }
        }
        probability_interp[n - 1] = _probability[_probability.size() - 1];
        sum_prob += probability_interp[n - 1];
    }

    std::vector<std::vector<double>> ranges;
    std::vector<double> current_range{0, 0, 0};
    std::vector<double> sorted_probabilities(probability_interp.begin(), probability_interp.end());
    double scaled_prob = probability * sum_prob;
    // Don't bother to store ranges with probability less than this
    const double min_prob = 0.005 * sum_prob;
    std::vector<int> perm(n);
    int num_values = 0;

    for (int i = 0; i < n; i++) perm[i] = i;
    revsort(&sorted_probabilities[0], &perm[0], n);

    double cumulative_prob= 0.;
    for (int i = 0; i < n; i++) {
        cumulative_prob += sorted_probabilities[i];
        num_values++;
        if (cumulative_prob > scaled_prob) {
            break;
        }
    }

    std::vector<int> included_values(perm.begin(), perm.begin() + num_values);
    std::sort(included_values.begin(), included_values.end());
    current_range[0] = _start_calAD + included_values[0] * resolution;
    current_range[2] = probability_interp[included_values[0]];
    for (int i = 1; i < num_values; i++) {
        if (included_values[i] - included_values[i-1] > 1) {
            current_range[1] = _start_calAD + included_values[i - 1] * resolution;
            current_range[2] /=  sum_prob;
            if (current_range[2] > min_prob) {
                ranges.push_back(current_range);
            }
            current_range[0] = _start_calAD + included_values[i] * resolution;
            current_range[2] = probability_interp[included_values[i]];
        } else {
            current_range[2] += probability_interp[included_values[i]];
        }
    }
    // last range
    if (current_range[2] > min_prob) {
        current_range[1] = _start_calAD + included_values[num_values - 1] * resolution;
        current_range[2] /= sum_prob;
        ranges.push_back(current_range);
    }
    return ranges;
}
