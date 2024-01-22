/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <iostream>
#include <fstream>
#include "DensityOutput.h"

/*
 * Creates an object suitable for printing density output in the format expected by the OxCal frontend.
 * Note that this class does not have enough information for printing a specific type of density, instead one of the
 * child classes PosteriorDensityOutput or PredictiveDensityOutput should be used instead.
 */
DensityOutput::DensityOutput(int index, double resolution)
    : _index(index), _resolution(resolution) {
    _output_var = "ocd[" + std::to_string(_index) + "]";
    _output_prefix = _output_var + ".posterior";
}

void DensityOutput::print() {
    std::string file_path = project_directory + project_name + ".js";
    std::ofstream output_file;
    output_file.open(file_path, std::ios_base::app);
    if (! output_file.is_open()) throw UnableToWriteToOutputFileException(file_path);

    for (const std::string& output_line : _get_output_lines()) output_file << output_line;
    output_file.close();
}

std::vector<std::string> DensityOutput::_get_output_lines() {
    std::vector<std::string> output_lines;
    int comment_index = 0;

    output_lines.push_back("if(!" + _output_var + "){" + _output_var + "={};}\n");
    output_lines.push_back(_variable_line("ref", carbondate_long_reference()));
    output_lines.push_back(_output_prefix + "={};\n");
    output_lines.push_back(_comment_line("Posterior ", comment_index));
    output_lines.push_back(_range_lines(comment_index));
    output_lines.push_back(_output_line("mean", _mean_calAD));
    output_lines.push_back(_output_line("sigma", _sigma_calAD));
    output_lines.push_back(_output_line("median", _median_calAD));
    output_lines.push_back(_output_line("probNorm", _prob_norm));
    output_lines.push_back(_output_line("start", _start_calAD));
    output_lines.push_back(_output_line("resolution", _resolution));
    output_lines.push_back(_output_line("prob", _probability));

    return output_lines;
}

std::string DensityOutput::_comment_line(const std::string& comment, int& comment_index) {
    std::string line;
    if (comment_index == 0) line = _output_prefix + ".comment=[];\n";
    line += _output_prefix + ".comment[" + std::to_string(comment_index++) + "]=\"" +  comment + "\";\n";
    return line;
}

std::string DensityOutput::_variable_line(const std::string& var_name, const std::string& var) {
    return _output_var + "." + var_name + "=\"" +  var + "\";\n";
}

std::string DensityOutput::_output_line(const std::string& var_name, double var) {
    return _output_prefix + "." + var_name + "=" + to_string(var, 6) + ";\n";
}

std::string DensityOutput::_output_line(const std::string& var_name, const std::vector<double>& var) {
    std::string output_line = _output_prefix + "." + var_name + "=[";
    for (int i = 0; i < var.size() - 1; i++) output_line += std::to_string(var[i]) + ", ";
    output_line += std::to_string(var[var.size() - 1]) + "];\n";
    return output_line;
}

void DensityOutput::_set_probability(const std::vector<double>& probability) {
    double prob_total = 0.;
    _probability.resize(probability.size() + 2, 0); // This ensures the beginning and final probability is zero
    for (int i = 0; i < probability.size(); i++) {
        _probability[i + 1] = probability[i];
        prob_total += probability[i];
        if (probability[i] > _prob_max) _prob_max = probability[i];
    }
    // Scale so that the maximum value of the density is 1.
    for (double & prob : _probability) prob /= _prob_max;
    _prob_norm = _prob_max / (prob_total * _resolution);

    _smoothed_probability.resize(_probability.size());
    int last = (int) _probability.size() - 1;
    _smoothed_probability[0] = (2. * _probability[0] + _probability[1]) / 3.;
    for (int i = 1; i < last; i++)
        _smoothed_probability[i] = (_probability[i - 1] + _probability[i] + _probability[i + 1]) / 3.;
    _smoothed_probability[last - 1] = (_probability[last - 1] + 2. * _probability[last]) / 3.;
}


std::string DensityOutput::_range_lines(int &comment_index) {
    // Must be implemented in child function if range lines are calculated.
    return {};
}
