#include <fstream>
#include <regex>
#include <set>
#include "read_data.h"
#include "csv_helpers.h"
#include "log.h"

#ifndef CALIBRATION_DATA_PREFIX
#define CALIBRATION_DATA_PREFIX "../oxcal/"
#endif

#ifndef DATA_PREFIX
#define DATA_PREFIX "../data/"
#endif

#ifndef OUTPUT_PREFIX
#define OUTPUT_PREFIX "../output/"
#endif

const std::set<std::string> modern_intcal_curves = {"intcal04.14c", "intcal09.14c", "intcal13.14c", "intcal20.14c"};
const std::set<std::string> old_intcal_curves = {"intcal98.14c"};
const std::set<std::string> custom_curves = {"HOBS2022.14c"};

std::string project_name;

void read_arguments(int argc, char* argv[]) {
    if (argc < 2)
        throw IncorrectArgumentsException();
    project_name = argv[1];
}

std::string oxcal_file_path() {
    return DATA_PREFIX + project_name + ".oxcal";
}

void read_calibration_curve(
        const std::string& calibration_curve,
        std::vector<double>& cc_cal_age,
        std::vector<double>& cc_c14_age,
        std::vector<double>& cc_c14_sig) {

    const std::string calibration_curve_path = CALIBRATION_DATA_PREFIX + calibration_curve;

    // Check we can read the file
    std::fstream file(calibration_curve_path, std::ios::in);
    if(!file.is_open()) throw UnableToReadCalibrationCurveException(calibration_curve_path);

    update_log_file("Reading calibration data from " + calibration_curve);
    if (modern_intcal_curves.count(calibration_curve) == 1) {
        cc_cal_age = get_csv_data_from_column(&file, 0, ',');
        cc_c14_age = get_csv_data_from_column(&file, 1, ',');
        cc_c14_sig = get_csv_data_from_column(&file, 2, ',');
    } else if (old_intcal_curves.count(calibration_curve) == 1) {
        cc_cal_age = get_csv_data_from_column(&file, 0, ' ');
        cc_c14_age = get_csv_data_from_column(&file, 3, ' ');
        cc_c14_sig = get_csv_data_from_column(&file, 4, ' ');
        for (double & cal_age : cc_cal_age) cal_age = 1950. - cal_age;
    } else if (custom_curves.count(calibration_curve) == 1) {
        cc_cal_age = get_csv_data_from_column(&file, 0, '\t');
        cc_c14_age = get_csv_data_from_column(&file, 3, '\t');
        cc_c14_sig = get_csv_data_from_column(&file, 4, '\t');
        for (double & cal_age : cc_cal_age) cal_age = 1950. - cal_age;
    }
}

// Takes a *.oxcal input file created by the OxCal software and reads it to determine the
// NP model data and output options.
// If NP model data is found it returns true, otherwise it returns false.
bool read_oxcal_data(
        std::vector<double>& c14_age,
        std::vector<double>& c14_sig,
        std::vector<double>& f14c_age,
        std::vector<double>& f14c_sig,
        std::string& model_name) {

    std::string line, age, sig, end_of_section = "};";
    std::regex np_model_regex(R"(NP_Model\(\s*["'](.*)["']\s*\))");
    std::regex named_r_date_regex(R"(R_Date\(\s*["'](.*)["']\s*,\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    std::regex unnamed_r_date_regex(R"(R_Date\(\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    std::regex named_r_f14c_regex(R"(R_F14C\(\s*["'](.*)["']\s*,\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    std::regex unnamed_r_f14c_regex(R"(R_F14C\(\s*([0-9\.]*)\s*,\s*([0-9\.]*))");
    bool np_model = false;
    std::smatch r_date_match;
    std::string filepath = oxcal_file_path();

    std::fstream file(filepath, std::ios::in);
    if(!file.is_open()) throw UnableToReadOxcalFileException(filepath);

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
        } else if (np_model && regex_search(line, r_date_match, named_r_f14c_regex)) {
            age = r_date_match[2];
            sig = r_date_match[3];
            f14c_age.push_back(std::strtod(age.c_str(), nullptr));
            f14c_sig.push_back(std::strtod(sig.c_str(), nullptr));
        } else if (np_model && regex_search(line, r_date_match, unnamed_r_f14c_regex)){
            age = r_date_match[1];
            sig = r_date_match[2];
            f14c_age.push_back(std::strtod(age.c_str(), nullptr));
            f14c_sig.push_back(std::strtod(sig.c_str(), nullptr));
        }
    }
    if (!c14_age.empty() && !f14c_age.empty())  throw InconsistentDateFormatsException();
    return np_model && (!c14_age.empty() || !f14c_age.empty());
}

int read_output_offset(const std::string& model_name) {
    std::string line, model_index;
    std::regex np_output_regex(R"(ocd\[([0-9]+)\].name\s*=\s*["'])" + model_name + R"(["'];)");
    std::smatch np_model_match;

    std::string output_file_path = OUTPUT_PREFIX + project_name + ".js";
    std::fstream file(output_file_path, std::ios::in);
    if(!file.is_open()) throw UnableToReadOutputFileException(output_file_path);

    while (getline(file, line)) {
        if (regex_search(line, np_model_match, np_output_regex)) {
            model_index = np_model_match[1];
            return std::stoi(model_index);
        }
    }
    throw UnableToDetermineOutputOffsetException(output_file_path, model_name);
}

// Reads in the options from the oxcal data file. If any of the options are not found in the file
// then the value will not be altered from the original value. Currently, it will populate the
// following option variables provided as arguments:
// * iterations: The number of iterations for the DPMM
// * resolution: The resolution used for outputting the predictive and posterior density
// * ranges: A vector of 3 values denoting whether to log the 68.3%, 95.4% and 99.7% ranges
// * quantile ranges: False to use HPD ranges, True to simply work out a single quantile range

void read_options(
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        bool &use_f14c,
        std::string &calibration_curve) {

    std::string line, option, value, end_of_section = "};";
    std::regex options_regex(R"(Options\(\s*\))");
    std::regex option_regex(R"((\w+)=['"]*([^'"]+)['"]*;)");
    bool options_block = false;
    std::smatch option_match;
    std::set<std::string> allowed_calibration_curves;
    allowed_calibration_curves.insert(modern_intcal_curves.begin(), modern_intcal_curves.end());
    allowed_calibration_curves.insert(old_intcal_curves.begin(), old_intcal_curves.end());
    allowed_calibration_curves.insert(custom_curves.begin(), custom_curves.end());

    std::string filepath = oxcal_file_path();
    std::fstream file(filepath, std::ios::in);
    if(!file.is_open()) throw UnableToReadOxcalFileException(filepath);

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
                ranges[0] = value == "TRUE";
            } else if (option == "SD2") {
                ranges[1] = value == "TRUE";
            } else if (option == "SD3") {
                ranges[2] = value == "TRUE";
            } else if (option == "Floruit") {
                quantile_ranges = value == "TRUE";
            } else if (option == "UseF14C") {
                use_f14c = value == "TRUE";
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
