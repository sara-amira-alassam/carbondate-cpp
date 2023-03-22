//
// Created by Sara Al-Assam on 13/03/2023.
//
#include <fstream>
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

//
bool read_oxcal_data(
        const std::string& file_prefix,
        std::vector<double>& c14_age,
        std::vector<double>& c14_sig,
        std::string& model_name) {

    std::string line, age, sig, end_of_section = "};";
    std::regex np_model_regex(R"(NP_Model\(\s*["'](.*)["']\s*\))");
    std::regex named_r_date_regex(R"(R_Date\(\s*["'](.*)["']\s*,\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    std::regex unnamed_r_date_regex(R"(R_Date\(\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    bool np_model = false;
    std::smatch r_date_match;
    std::fstream file("../data/" + file_prefix + ".oxcal", std::ios::in);

    if(!file.is_open()) throw std::runtime_error("Could not open file");

    while (getline(file, line)) {
        if (regex_search(line, r_date_match, np_model_regex)) {
            np_model = true;
            model_name = r_date_match[1];
        } else if (np_model && line.find("};") != std::string::npos) {
            break;
        } else if (np_model && regex_search(line, r_date_match, named_r_date_regex)) {
            age = r_date_match[2];
            sig = r_date_match[3];
            c14_age.push_back(std::strtod(age.c_str(), nullptr));
            c14_sig.push_back(std::strtod(sig.c_str(), nullptr));
        } else if (np_model && regex_search(line, r_date_match, unnamed_r_date_regex)){
            age = r_date_match[1];
            sig = r_date_match[2];
            c14_age.push_back(std::strtod(age.c_str(), nullptr));
            c14_sig.push_back(std::strtod(sig.c_str(), nullptr));
        }
    }
    return np_model && !c14_age.empty();
}
