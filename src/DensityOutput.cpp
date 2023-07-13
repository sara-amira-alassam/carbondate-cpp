#include <iostream>
#include <fstream>
#include "DensityOutput.h"

DensityOutput::DensityOutput(int index, double resolution)
    : _index(index), _resolution(resolution) {
    _output_var = "ocd[" + std::to_string(_index) + "]";
    _output_prefix = _output_var + ".posterior";
}

void DensityOutput::print() {
    std::string file_path = "../output/" + project_name + ".js";
    std::ofstream output_file;
    output_file.open(file_path, std::ios_base::app);
    if (! output_file.is_open()) throw UnableToWriteToOutputFileException(file_path);

    for (const std::string& output_line : get_output_lines()) output_file << output_line;
    output_file.close();
}

std::vector<std::string> DensityOutput::get_output_lines() {
    std::vector<std::string> output_lines;
    int comment_index = 0;

    output_lines.push_back("if(!" + _output_var + "){" + _output_var + "={};}\n");
    output_lines.push_back(variable_line("ref", carbondate_long_reference()));
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
    if (comment_index == 0) line = _output_prefix + ".comment=[];\n";
    line += _output_prefix + ".comment[" + std::to_string(comment_index++) + "]=\"" +  comment + "\";\n";
    return line;
}

std::string DensityOutput::variable_line(const std::string& var_name, const std::string& var) {
    return _output_var + "." + var_name + "=\"" +  var + "\";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, double var) {
    return _output_prefix + "." + var_name + "=" + to_string(var, 6) + ";\n";
}

std::string DensityOutput::output_line(const std::string& var_name, const std::vector<double>& var) {
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
    // Scale so that the maximum value of the density is 1.
    for (int i = 0; i < probability.size(); i++) _probability[i] /= _prob_max;
    _prob_norm = _prob_max / (prob_total * _resolution);
}


std::string DensityOutput::range_lines(int &comment_index) {
    // Must be implemented in child function if range lines are calculated.
    return {};
}
