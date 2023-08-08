#include <fstream>
#include "plain_text.h"

#ifdef OXCAL_RELEASE
#define TEXT_PREFIX ""
#else
#define TEXT_PREFIX "../output/"
#endif

std::string text_file_path() {
    return TEXT_PREFIX + project_name + ".txt";
}

void initialize_text_file() {
    std::string filepath = text_file_path();
#ifdef OXCAL_RELEASE  // We only expect the file to exist already if running on OxCal
    std::ofstream file(filepath, std::fstream::in);
    if (!file.is_open()) throw UnableToFindTextFileException(filepath);

    file.close();
#endif
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