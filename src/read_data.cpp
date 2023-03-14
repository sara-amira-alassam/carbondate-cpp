//
// Created by Sara Al-Assam on 13/03/2023.
//
#include <fstream>
#include <sstream>
#include <regex>
#include "read_data.h"
#include "csv_helpers.h"

void read_calibration_curve(
        const std::string& calibration_curve_path,
        std::vector<double>& cc_cal_age,
        std::vector<double>& cc_c14_age,
        std::vector<double>& cc_c14_sig) {

    cc_cal_age = get_csv_data_from_column(calibration_curve_path, 0);
    cc_c14_age = get_csv_data_from_column(calibration_curve_path, 1);
    cc_c14_sig = get_csv_data_from_column(calibration_curve_path, 2);

}

void read_oxcal_data(
        const std::string& file_prefix,
        std::vector<std::string>& c14_name,
        std::vector<double>& c14_age,
        std::vector<double>& c14_sig) {

    double val;
    std::string line, name, age, sig;
    std::regex r_date_regex(R"(R_Date\(\s*["'](.*)["']\s*,\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    std::smatch r_date_match;
    bool is_r_date;

    std::fstream file("../data/" + file_prefix + ".oxcal", std::ios::in);

    if(!file.is_open()) throw std::runtime_error("Could not open file");

    while (getline(file, line)) {
        is_r_date = regex_search(line, r_date_match, r_date_regex);
        if (is_r_date && r_date_match.size() == 4) {
            name = r_date_match[1];
            age = r_date_match[2];
            sig = r_date_match[3];
            c14_name.push_back(name);
            c14_age.push_back(std::strtod(age.c_str(), nullptr));
            c14_sig.push_back(std::strtod(sig.c_str(), nullptr));
        }
    }
}
