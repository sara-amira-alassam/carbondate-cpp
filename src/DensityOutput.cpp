//
// Created by Sara Al-Assam on 21/02/2023.
//
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include "helpers.h"
#include "DensityOutput.h"
#include <iostream>

DensityOutput::DensityOutput(
        std::string output_var,
        int index,
        const std::string& output_name,
        double date,
        double error,
        const std::string& name) {
    _output_var = std::move(output_var) + "[" + std::to_string(index) + "]";
    _output_prefix = _output_var + "." + output_name;
    _date = date;
    _error = error;
    _name = name;
}

void DensityOutput::write_to_file(int resolution, const std::string &file_prefix) {
    std::ofstream output_file;
    output_file.open("../output/" + file_prefix + ".js", std::ios_base::app);

    for (const std::string& output_line : get_output_lines(resolution)) {
        output_file << output_line;
    }

    output_file.close();
}

void DensityOutput::print(int resolution) {
    for (const std::string& output_line : get_output_lines(resolution)) {
        std::cout << output_line;
    }
}

std::vector<std::string> DensityOutput::get_output_lines(int resolution) {
    std::vector<std::string> output_lines;
    int comment_index = 0, no_comments = -1;
    std::string param = to_string(_date) + "," + to_string(_error);
    calculate_probability_smoothed(resolution);

    output_lines.push_back("if(!" + _output_var + "){" + _output_var + "={};}\n");
    output_lines.push_back(variable_line("ref", "OxCal v4.4.4 Bronk Ramsey (2021); r:5"));
    output_lines.push_back(variable_line("name", _name));
    output_lines.push_back(variable_line("op", "R_Date"));
    output_lines.push_back(variable_line("param", param));
    output_lines.push_back(variable_line("level", 0));
    output_lines.push_back(variable_line("date", _date));
    output_lines.push_back(variable_line("error", _error));
    output_lines.push_back(_output_prefix + "={};\n");
    output_lines.push_back(comment_line(_name + " R_Date(" + param + ")", comment_index));
    output_lines.push_back(variable_line("type", "date"));
    output_lines.push_back(variable_line("calib", 0));
    output_lines.push_back(range_lines(1, 0.683, resolution, comment_index));
    output_lines.push_back(range_lines(2, 0.954, resolution, comment_index));
    output_lines.push_back(range_lines(3, 0.997, resolution, no_comments));
    output_lines.push_back(output_line("mean", mean_calAD));
    output_lines.push_back(output_line("sigma", sigma));
    output_lines.push_back(output_line("median", median_calAD));
    output_lines.push_back(output_line("probNorm", _prob_norm_smoothed));
    output_lines.push_back(output_line("start", _start_calAD_smoothed));
    output_lines.push_back(output_line("resolution", resolution));
    output_lines.push_back(output_line("prob", _prob_smoothed));

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

std::string DensityOutput::variable_line(const std::string& var_name, int var) {
    return _output_var + "." + var_name + "=" +  std::to_string(var) + ";\n";
}

std::string DensityOutput::variable_line(const std::string& var_name, double var) {
    return _output_var + "." + var_name + "=" +  to_string(var) + ";\n";
}

std::string DensityOutput::variable_line(const std::string& var_name, const std::string& var) {
    return _output_var + "." + var_name + "=\"" +  var + "\";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, int var) {
    return _output_prefix + "." + var_name + "=" +  std::to_string(var) + ";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, double var) {
    return _output_prefix + "." + var_name + "=" + to_string(var) + ";\n";
}

std::string DensityOutput::output_line(
        const std::string& var_name, const std::vector<double>& var) {
    std::string output_line = _output_prefix + "." + var_name + "=[";
    for (int i = 0; i < var.size() - 1; i++) output_line += std::to_string(var[i]) + ", ";
    output_line += std::to_string(var[var.size() - 1]) + "];\n";
    return output_line;
}

std::string DensityOutput::range_lines(
        int range_index, double probability, int resolution, int& comment_index) {

    std::vector<std::vector<double>> ranges = get_ranges_by_bisection(probability, resolution);
    std::string range_string = "range[" + std::to_string(range_index) + "]";
    std::string range_lines;
    if (range_index == 1) {
        range_lines = _output_prefix + ".range=[];\n";
    }
    range_lines += _output_prefix + "." + range_string + "=[];\n";
    for (int i = 0; i < ranges.size(); i++) {
        range_lines += output_line(range_string + "[" + std::to_string(i) + "]", ranges[i]);
    }
    if (comment_index >= 0) {
        range_lines += comment_line(
                "  " + to_percent_string(probability) + " probability", comment_index);
        for (auto & range : ranges) {
            std::string comment = "    " + std::to_string(int (round(range[0]))) + "AD";
            comment += " (" + to_percent_string(range[2]) + ") ";
            comment += std::to_string(int (round(range[1]))) + "AD";
            range_lines += comment_line(comment, comment_index);
        }
    }
    return range_lines;
}

std::string DensityOutput::to_string(double var) {
    std::string string_var;
    char temp_string[10] = "";
    snprintf(temp_string, 10, "%.6g", var);
    string_var = temp_string;
    return string_var;
}

std::string DensityOutput::to_percent_string(double fraction) {
    std::string percent;
    char temp_string[8] = "";
    snprintf(temp_string, 8, "%4.1f%%", fraction * 100);
    percent = temp_string;
    return percent;
}

std::vector<std::vector<double>> DensityOutput::as_columns(int resolution) {
    calculate_probability_smoothed(resolution);
    std::vector<std::vector<double>> output(2, std::vector<double>(_prob_smoothed.size()));
    for (int i = 0; i < _prob_smoothed.size(); i++) {
        output[0][i] = 1950.5 - _start_calAD_smoothed - i * resolution;
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
    _prob_norm_smoothed = max_density / resolution;
}

void DensityOutput::set_yearwise_probability(std::vector<double> probability) {
    double prob_total = 0.;
    _prob_yearwise.resize(probability.size());
    _prob_max = 0.;
    // We also need to clear variables related to the smoothed probability
    _prob_smoothed.clear();
    _resolution_smoothed = 0;
    for (int i = 0; i < probability.size(); i++) {
        _prob_yearwise[i] = probability[i];
        prob_total += probability[i];
        if (probability[i] > _prob_max) _prob_max = probability[i];
    }
    for (int i = 0; i < probability.size(); i++) _prob_yearwise[i] /= prob_total;
}

// Returns the area under the yearwise probability curve if we ignore all values below the cut-off
// Also populates the vector ranges, where each entry contains
// [start_calAD, end_calAD, probability within this range]
double DensityOutput::find_probability_and_ranges_for_cut_off_smoothed(
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
            ranges.push_back(
                    std::vector<double> {x_intercept_1, x_intercept_2, 100. * range_probability});
            total_probability += range_probability;
            range_probability = 0;
        } else if (cut_off <= y1 && cut_off <= y2) {
            range_probability += (y1 + y2) * res / 2.;
        }
    }
    return total_probability;
}

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges_by_bisection(
        double probability, int resolution) {
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
        current_probability = find_probability_and_ranges_for_cut_off_smoothed(p, ranges);
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

/////////////////// NOTE: none of these functions are currently used but are left in for testing
/////////////////// They should be removed when the final version is decided on.
// Returns the area under the yearwise probability curve if we ignore all values below the cut-off
// Also populates the vector ranges, where each entry contains
// [start_calAD, end_calAD, probability within this range]
/*
double DensityOutput::find_probability_and_ranges_for_cut_off(
        double cut_off, std::vector<std::vector<double>> &ranges) {
    ranges.clear();
    double y1, y2, dx, x_intercept_1, x_intercept_2;
    double range_probability = 0, total_probability = 0;
    for (int i = 0; i < _prob_yearwise.size() - 1; i++) {
        y1 = _prob_yearwise[i];
        y2 = _prob_yearwise[i + 1];
        if (y1 <= cut_off && cut_off <= y2) {
            dx = (cut_off - y1) / (y2 - y1);
            x_intercept_1 = start_calAD + i + dx;
            range_probability = (cut_off + y2) * (1 - dx) / 2.;
        } else if (y2 <= cut_off && cut_off <= y1) {
            dx = (cut_off - y1) / (y2 - y1);
            x_intercept_2 = start_calAD + i + dx;
            range_probability += (cut_off + y2) * dx / 2.;
            ranges.push_back(std::vector<double> {x_intercept_1, x_intercept_2, range_probability});
            total_probability += range_probability;
            range_probability = 0;
        } else if (cut_off <= y1 && cut_off <= y2) {
            range_probability += (y1 + y2) / 2.;
        }
    }
    return total_probability;
}

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges_by_bisection(double probability) {
    std::vector<std::vector<double>> ranges;

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

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the ressult as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> DensityOutput::get_ranges(double probability) {

    std::vector<std::vector<double>> ranges;
    std::vector<double> current_range{0, 0, 0};
    std::vector<double> sorted_probabilities(_prob_yearwise.begin(), _prob_yearwise.end());
    double scaled_prob = probability;
    const double min_prob = 0.005; // Don't bother to store ranges with probability less than this
    int n = (int) sorted_probabilities.size();
    std::vector<int> perm(n);
    int num_values = 0;

    for (int i = 0; i < n; i++) perm[i] = i;
    revsort(&sorted_probabilities[0], &perm[0], n);

    double cumulative_prob= 0.;
    for (int i = 0; i < n; i++) {
        cumulative_prob += sorted_probabilities[i];
        num_values++;
        if (abs(cumulative_prob - scaled_prob) < 1e-3 || cumulative_prob > scaled_prob) {
            break;
        }
    }

    std::vector<int> included_values(perm.begin(), perm.begin() + num_values);
    std::sort(included_values.begin(), included_values.end());
    current_range[0] = start_calAD + included_values[0];
    current_range[2] = _prob_yearwise[included_values[0]];
    for (int i = 1; i < num_values; i++) {
        if (included_values[i] - included_values[i-1] > 1) {
            current_range[1] = start_calAD + included_values[i - 1];
            if (current_range[2] > min_prob) {
                ranges.push_back(current_range);
            }
            current_range[0] = start_calAD + included_values[i];
            current_range[2] = _prob_yearwise[included_values[i]];
        } else {
            current_range[2] += _prob_yearwise[included_values[i]];
        }
    }
    // last range
    if (current_range[2] > min_prob) {
        current_range[1] = start_calAD + included_values[num_values - 1];
        ranges.push_back(current_range);
    }
    return ranges;
}
*/
