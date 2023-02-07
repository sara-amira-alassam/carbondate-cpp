#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

std::vector<double> get_csv_data_from_column(const std::string& filename, int column_index) {

    std::vector<double> data_column;
    double val;
    std::string line, word;
    char *ptr;

    std::fstream file(filename, std::ios::in);

    if(!file.is_open()) throw std::runtime_error("Could not open file");

    while (getline(file, line)) {
        std::stringstream str(line);
       int ind = 0;
        while (getline(str, word, ',')) {
            if (ind++ == column_index) {
                val = std::strtod(word.c_str(), &ptr);
                if (*ptr == '\0') data_column.push_back(val);
                break;
            }
        }
    }
    return data_column;
}

void write_columns_to_csv(
        const std::string& filename,
        std::vector<std::string> headers,
        std::vector<std::vector<double>> data) {
    // TODO: Should check all the columns are the same size
    std::ofstream file;
    file.open(filename);

    char header_line[256];
    for (int j = 0; j < data.size(); j++) {
        if (j == 0) {
            sprintf(header_line, "%s", headers[j].c_str());
        } else {
            sprintf(header_line, "%s%s", header_line, headers[j].c_str());
        }
        if (j < data.size() - 1) sprintf(header_line, "%s, ", header_line);
    }
    sprintf(header_line, "%s\n", header_line);
    file << header_line;

    for (int i = 0; i < data[0].size(); i++) {
        char line[256];
        for (int j = 0; j < data.size(); j++) {
            if (j == 0) {
                sprintf(line, "%.10e", data[j][i]);
            } else {
                sprintf(line, "%s%.10e", line, data[j][i]);
            }
            if (j < data.size() - 1) sprintf(line, "%s, ", line);
        }
        sprintf(line, "%s\n", line);
        file << line;
    }
}
