//
// Created by Sara Al-Assam on 21/02/2023.
//
#include <utility>
#include <cmath>
#include "DensityOutput.h"
#include <iostream>

DensityOutput::DensityOutput(std::string output_var, int index, std::string output_name) {
    _output_var = std::move(output_var) + "[" + std::to_string(index) + "]";
    _index = index;
    _output_name = std::move(output_name);
}

std::string DensityOutput::output_prefix() {
    if (_output_prefix.empty()) {
        _output_prefix = _output_var + "." + _output_name;
    }
    return _output_prefix;
}

void DensityOutput::print(int resolution) {
    std::cout << "if!(" + _output_var + "){" + _output_var + "={};}\n";
    std::cout << output_line("mean", mean_calAD);
    std::cout << output_line("sigma", sigma);
    std::cout << output_line("median", median_calAD);
    get_ranges(0.683, resolution);
    get_ranges(0.954, resolution);
    get_ranges(0.997, resolution);
    std::cout << output_line("start", _start_calAD_smoothed);
    std::cout << output_line("resolution", resolution);
    std::cout << output_line("prob", _prob_smoothed);
    std::cout << output_line("probNorm", _prob_norm_smoothed);
}

std::string DensityOutput::output_line(const std::string& var_name, int var) {
    return output_prefix() + "." + var_name + "=" +  std::to_string(var) + ";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, double var) {
    return output_prefix() + "." + var_name + "=" +  std::to_string(var) + ";\n";
}

std::string DensityOutput::output_line(
        const std::string& var_name, const std::vector<double>& var) {
    std::string output_line = output_prefix() + "." + var_name + "=[";
    for (int i = 0; i < var.size() - 1; i++) output_line += std::to_string(var[i]) + ", ";
    output_line += std::to_string(var[var.size() - 1]) + "];\n";
    return output_line;
}

std::vector<std::vector<double>> DensityOutput::as_columns(int resolution) {
    calculate_probability_smoothed(resolution);
    std::vector<std::vector<double>> output(2, std::vector<double>(_prob_smoothed.size()));
    for (int i = 0; i < _prob_smoothed.size(); i++) {
        output[0][i] = 1950 - _start_calAD_smoothed - i * resolution;
        output[1][i] = _prob_smoothed[i] * _prob_norm_smoothed;
    }
    return output;
}

void DensityOutput::calculate_probability_smoothed(int resolution) {
    if (!_prob_smoothed.empty() && _resolution_smoothed == resolution) {
        // This means its already been previously calculated and set
        return;
    }
    int break_num, num_breaks = (int) ceil((double) _prob_yearwise.size() / resolution);
    _prob_smoothed.resize(num_breaks, 0);
    _resolution_smoothed = resolution;
    _start_calAD_smoothed = start_calAD + (resolution - 1.) / 2;
    double max_density = 0.;
    for (int j = 0; j < _prob_yearwise.size(); j++) {
        break_num = j / resolution;
        _prob_smoothed[break_num] += _prob_yearwise[j];
        if (_prob_smoothed[break_num] > max_density) max_density = _prob_smoothed[break_num];
    }
    for (int b = 0; b < num_breaks; b++) _prob_smoothed[b] /= max_density;
    _prob_norm_smoothed = max_density / (_prob_total * resolution);
}

void DensityOutput::set_yearwise_probability(std::vector<double> probability) {
    _prob_yearwise.resize(probability.size());
    _prob_total = 0.;
    _prob_max = 0.;
    // We also need to clear variables related to the smoothed probability
    _prob_smoothed.clear();
    _resolution_smoothed = 0;
    for (int i = 0; i < probability.size(); i++) {
        _prob_yearwise[i] = probability[i];
        _prob_total += probability[i];
        if (probability[i] > _prob_max) _prob_max = probability[i];
    }
}

// Returns the area under the yearwise probability curve if we ignore all values below the cut-off
// Also populates the vector ranges, where each entry contains
// [start_calAD, end_calAD, probability within this range]
double DensityOutput::find_probability_and_ranges_for_cut_off(
        double cut_off, std::vector<std::vector<double>> &ranges) {
    ranges.clear();
    double y1, y2, dx, x_intercept_1, x_intercept_2, res = _resolution_smoothed;
    double range_probability = 0, total_probability = 0;
    for (int i = 0; i < _prob_smoothed.size() - 1; i++) {
        y1 = _prob_smoothed[i];
        y2 = _prob_smoothed[i + 1];
        if (y1 <= cut_off && cut_off <= y2) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_1 = _start_calAD_smoothed + i * res + dx;
            range_probability = (cut_off + y2) * (res - dx) / 2.;
        } else if (y2 <= cut_off && cut_off <= y1) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_2 = _start_calAD_smoothed + i * res + dx;
            range_probability += (cut_off + y2) * dx / 2.;
            range_probability *= _prob_norm_smoothed;
            ranges.push_back(std::vector<double> {x_intercept_1, x_intercept_2, range_probability});
            total_probability += range_probability;
            range_probability = 0;
        } else if (cut_off <= y1 && cut_off <= y2) {
            range_probability += (y1 + y2) * res / 2.;
        }
    }
    return total_probability;
}

// Finds the calendar ange ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges(double probability, int resolution) {
    std::vector<std::vector<double>> ranges;

    calculate_probability_smoothed(resolution);

    // Use bisection method to find the closest probability cut-off to give the desired probability
    // a and b are the upper and lower points of the section - we know the smoothed probability has a max of 1
    double a = 0., b = 1.;
    double p; // p is the midpoint between a and b
    double current_probability = -1.;
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
