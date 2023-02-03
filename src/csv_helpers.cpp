#include <fstream>
#include <vector>
#include <sstream>

std::vector<double> get_csv_data_from_column(const std::string& filename, unsigned column_index) {

    std::vector<double> data_column;
    double val;
    std::string line, word;
    char *ptr;

    std::fstream file(filename, std::ios::in);

    if(!file.is_open()) throw std::runtime_error("Could not open file");

    while (getline(file, line)) {
        std::stringstream str(line);
        unsigned ind = 0;
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
