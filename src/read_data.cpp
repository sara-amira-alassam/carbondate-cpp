#include <fstream>
#include <regex>
#include <set>
#include "read_data.h"
#include "csv_helpers.h"



void read_calibration_curve(
        const std::string& calibration_curve_path,
        std::vector<double>& cc_cal_age,
        std::vector<double>& cc_c14_age,
        std::vector<double>& cc_c14_sig) {

    printf("Reading calibration data from %s\n", calibration_curve_path.c_str());
    cc_cal_age = get_csv_data_from_column(calibration_curve_path, 0);
    cc_c14_age = get_csv_data_from_column(calibration_curve_path, 1);
    cc_c14_sig = get_csv_data_from_column(calibration_curve_path, 2);

}

// Takes a *.oxcal input file created by the OxCal software and reads it to determine the
// NP model data and output options.
// If NP model data is found it returns true, otherwise it returns false.
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
        } else if (np_model && line.find(end_of_section) != std::string::npos) {
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

int read_output_offset(const std::string& file_prefix, const std::string& model_name) {
    std::string line, model_index;
    std::regex np_output_regex(R"(ocd\[([0-9]+)\].name\s*=\s*["'])" + model_name + R"(["'];)");
    std::smatch np_model_match;

    std::fstream file("../output/" + file_prefix + ".js", std::ios::in);
    while (getline(file, line)) {
        if (regex_search(line, np_model_match, np_output_regex)) {
            model_index = np_model_match[1];
            return std::stoi(model_index);
        }
    }
    throw std::runtime_error("Could not find NP model " + model_name + " in output file");
}

// Reads in the options from the oxcal data file. If any of the options are not found in the file
// then the value will not be altered from the original value. Currently, it will populate the
// following option variables provided as arguments:
// * iterations: The number of iterations for the DPMM
// * resolution: The resolution used for outputting the predictive and posterior density
// * ranges: A vector of 3 values denoting whether to calculate the 68.3%, 95.4% and 99.7% ranges
void read_options(
        const std::string &file_prefix,
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        bool &intercept_ranges,
        std::string &calibration_curve) {

    std::string line, option, value, end_of_section = "};";
    std::regex options_regex(R"(Options\(\s*\))");
    std::regex option_regex(R"((\w+)=['"]*([^'"]+)['"]*;)");
    bool options_block = false;
    std::smatch option_match;
    std::fstream file("../data/" + file_prefix + ".oxcal", std::ios::in);
    static const std::set<std::string> allowed_calibration_curves{
        "intcal04.14c", "intcal09.14c", "intcal13.14c", "intcal20.14c"};

    if(!file.is_open()) throw std::runtime_error("Could not open file");

    while (getline(file, line)) {
        if (regex_search(line, option_match, options_regex)) {
            options_block = true;
        } else if (options_block && line.find(end_of_section) != std::string::npos) {
            break;
        } else if (options_block && regex_search(line, option_match, option_regex)) {
            option = option_match[1];
            value = option_match[2];
            if (option == "Resolution") {
                resolution = std::stod(value);
            } else if (option == "kIterations") {
                iterations = std::stoi(value) * 1000;
            } else if (option == "SD1") {
                ranges[0] = value != "FALSE";
            } else if (option == "SD2") {
                ranges[1] = value == "TRUE";
            } else if (option == "SD3") {
                ranges[2] = value == "TRUE";
            } else if (option == "Floruit") {
                quantile_ranges = value == "TRUE";
            } else if (option == "Intercept") {
                intercept_ranges = value == "TRUE";
            } else if (option == "Curve") {
                if (allowed_calibration_curves.count(value) == 0) {
                    printf(
                        "The calibration curve %s is not recognised. Default %s is being used.\n",
                        value.c_str(), calibration_curve.c_str());
                } else {
                    calibration_curve = value;
                }
            }
        }
    }
}
