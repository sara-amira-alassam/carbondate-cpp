#include <fstream>
#include "plain_text.h"

std::string text_file_path() {
    return project_directory + project_name + ".txt";
}

void initialize_text_file() {
    std::string filepath = text_file_path();
}

void update_text_file(const std::string& label, const std::vector<double>& values) {
    std::string filepath = text_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToTextFileException(filepath);

    std::string text_line = "@" + label;

    for (double value : values) text_line += "\t" + to_string(value, 4);

    file << text_line << std::endl;
    file.close();
}
