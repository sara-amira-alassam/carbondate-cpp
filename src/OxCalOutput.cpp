//
// Created by Sara Al-Assam on 15/03/2023.
//
#include <fstream>
#include <utility>
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

OxCalOutput::OxCalOutput(std::string file_prefix): _file_prefix {std::move(file_prefix)} {
    initialise_file();
}
