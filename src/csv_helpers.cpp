#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

std::vector<double> get_csv_data_from_column(std::fstream* file, int column_index, const char separator) {

    std::vector<double> data_column;
    double val;
    std::string line, word;
    char *ptr;

    file->clear();
    file->seekg(0, std::ios_base::beg);

    while (getline(*file, line)) {
        if (line[0] == '!' or line[0] == '#')  // comment line
            continue;
        std::stringstream str(line);
        int ind = 0;
        while (getline(str, word, separator)) {
            if (separator == ' ' && word.length() == 0)  // accounting for multiple whitespace is space separator
                continue;
            if (ind++ == column_index) {
                val = std::strtod(word.c_str(), &ptr);
                if (*ptr == '\0') data_column.push_back(val);
                break;
            }
        }
    }
    return data_column;
}
