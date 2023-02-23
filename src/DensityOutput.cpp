//
// Created by Sara Al-Assam on 21/02/2023.
//
#include <utility>
#include "DensityOutput.h"
#include <iostream>

DensityOutput::DensityOutput(std::string output_var, int index, std::string output_name) {
    _output_var = std::move(output_var);
    _index = index;
    _output_name = std::move(output_name);
}

std::string DensityOutput::output_prefix() {
    if (_output_prefix.empty()) {
        _output_prefix = _output_var + "[" + std::to_string(_index) + "]." + _output_name;
    }
    return _output_prefix;
}

void DensityOutput::print() {
    std::cout << output_line("mean", mean_calAD);
    std::cout << output_line("sigma", sigma);
    std::cout << output_line("median", median_calAD);
    std::cout << output_line("probNorm", prob_norm);
    std::cout << output_line("start", start_calAD);
    std::cout << output_line("resolution", resolution);
    std::cout << output_line("prob", prob);
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

std::vector<std::vector<double>> DensityOutput::as_columns() {
    std::vector<std::vector<double>> output(2, std::vector<double>(prob.size()));
    for (int i = 0; i < prob.size(); i++) {
        output[0][i] = start_calAD + i * resolution;
        output[1][i] = prob[i] * prob_norm;
    }
    return output;
}
