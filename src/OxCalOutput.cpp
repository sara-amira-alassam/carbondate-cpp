//
// Created by Sara Al-Assam on 15/03/2023.
//
#include <fstream>
#include "helpers.h"
#include "OxCalOutput.h"

void OxCalOutput::initialise_file() {
    // Here we are manually copying the file that we expect to be returned from OxCal
    // This can be removed eventually as the file should be produced by the OxCal program
    // This is just a useful shortcut for testing
    std::ifstream input_file;
    std::ofstream output_file;
    std::string line;
    input_file.open("../output/sum.js");
    output_file.open("../output/" + _file_prefix + ".js", std::ios_base::trunc);

    if (input_file && output_file) {
        while (getline(input_file, line)) {
            output_file << line << "\n";
        }
    }
    input_file.close();
    output_file.close();
}

void OxCalOutput::print_model() {
    std::vector<std::string> output_lines;
    std::string np_line;
    // This is a hack for now until we get it working properly
    np_line = R"(model.element[1]= {op: "NP_Model", type: "date", name: "", )";
    np_line += "pos: 1, timepos: 1, ";
    np_line += "lower: 650.5, upper: 1190.5};\n";

    output_lines.emplace_back(np_line);
    append_to_file(output_lines);
}

void OxCalOutput::print_predictive_density() {

}

void OxCalOutput::print_posteriors() {
    for (int i = 0; i < _posteriors.size(); i++) {
        _posteriors[i].write_to_file(
                _resolution, _file_prefix, "ocd[" + std::to_string(i + 2) + "]", "posterior");
    }
}

void OxCalOutput::append_to_file(const std::vector<std::string>& output_lines) {
    std::ofstream output_file;

    output_file.open("../output/" + _file_prefix + ".js", std::ios_base::app);
    for (const std::string& output_line : output_lines) output_file << output_line;
    output_file.close();
}

OxCalOutput::OxCalOutput(int n_obs, int resolution, const std::string &file_prefix) {
    _resolution = resolution;
    _file_prefix = file_prefix;
    _posteriors.reserve(n_obs);
}

void OxCalOutput::set_posterior(int ident, const DensityOutput& posterior) {
    if (ident >= _posteriors.size()) {
        _posteriors.push_back(posterior);
    } else {
        _posteriors[ident] = posterior;
    }
}

void OxCalOutput::print() {
    initialise_file();
    print_model();
    print_predictive_density();
    print_posteriors();
}
