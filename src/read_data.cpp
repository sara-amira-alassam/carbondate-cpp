/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <fstream>
#include <sstream>
#include <regex>
#include <set>
#include "read_data.h"
#include "csv_helpers.h"
#include "log.h"

std::string calling_directory;
std::string project_name;
std::string project_directory;
std::string oxcal_version;

std::string get_path(const std::string& full_path) {

    size_t last_slash_pos = full_path.rfind('/');
    std::string path = full_path.substr(0, last_slash_pos + 1);

    return path;
}

/* Must be called before any other calls - this sets the global variable `project_name` which is used to write to
 * the various output files during execution, and to know which input file to read */
void read_arguments(int argc, char* argv[]) {
    if (argc < 2)
        throw IncorrectArgumentsException();
    calling_directory = get_path(argv[0]);
    project_directory = get_path(argv[1]);

    project_name = argv[1];
    std::regex project_name_regex(R"(/([^/]+)\.oxcal$)");
    std::smatch matches;

    if (regex_search(project_name, matches, project_name_regex)) {
        if (matches.size() > 1) project_name = matches[1].str();
    } else {
        throw IncorrectArgumentsException();
    }
}

std::string oxcal_file_path() {
    return project_directory + project_name + ".oxcal";
}

/*
 * Reads a calibration curve from a file, given the name of the calibration curve.
 * Args:
 *  - calibration_curve_name: One of the allowed calibration curve names
 *  - cc_cal_age: An empty vector to hold the calibration curve calendar ages
 *  - cc_c14_age: An empty vector to hold the calibration curve radiocarbon ages
 *  - cc_c14_sig: An empty vector to hold the radiocarbon age uncertainties
 */
void read_calibration_curve(
        const std::string& calibration_curve_name,
        std::vector<double>& cc_cal_age,
        std::vector<double>& cc_c14_age,
        std::vector<double>& cc_c14_sig) {

#ifdef OXCAL_RELEASE
    const std::string calibration_curve_path = calling_directory + calibration_curve_name;
#else
    const std::string calibration_curve_path = "../curves/" + calibration_curve_name;
#endif
    // Check we can read the file
    std::fstream file(calibration_curve_path, std::ios::in);
    if(!file.is_open()) throw UnableToReadCalibrationCurveException(calibration_curve_path);

    update_log_file("Reading calibration data from " + calibration_curve_name);
    if (calibration_curve_name == "intcal98.14c") {
        cc_cal_age = get_csv_data_from_column(&file, 0, ' ');
        cc_c14_age = get_csv_data_from_column(&file, 3, ' ');
        cc_c14_sig = get_csv_data_from_column(&file, 4, ' ');
        for (double & cal_age : cc_cal_age) cal_age = 1950. - cal_age;
    } else {
        cc_cal_age = get_csv_data_from_column(&file, 0, ',');
        cc_c14_age = get_csv_data_from_column(&file, 1, ',');
        cc_c14_sig = get_csv_data_from_column(&file, 2, ',');
    }
}

/*
 * Takes a *.oxcal input file created by the OxCal software and reads it to determine the NP model data.
 * If NP model data is found it returns true, otherwise it returns false.
 * Args:
 * - date_name: An empty string vector to hold the names of the data points
 * - c14_age: An empty vector to hold the radiocarbon ages of the data points. Will be empty if F14C ages are specified
 * - c14_sig: An empty vector to hold the radiocarbon age error. Will be empty if F14C ages are specified
 * - f14c_age: An empty vector to hold the radiocarbon ages of the data points. Will be empty if c14 ages are specified
 * - f14c_sig: An empty vector to hold the radiocarbon age error. Will be empty if c14 ages are specified
 * - model_name: The name given by the user to the NP model.
 */
bool read_oxcal_data(
        std::vector<std::string>& date_name,
        std::vector<double>& c14_age,
        std::vector<double>& c14_sig,
        std::vector<double>& f14c_age,
        std::vector<double>& f14c_sig,
        std::string& model_name) {

    std::string line, name, age, sig, end_of_section = "};";
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
            name = r_date_match[1];
            age = r_date_match[2];
            sig = r_date_match[3];
            date_name.push_back(name);
            c14_age.push_back(std::strtod(age.c_str(), nullptr));
            c14_sig.push_back(std::strtod(sig.c_str(), nullptr));
        } else if (np_model && regex_search(line, r_date_match, unnamed_r_date_regex)){
            age = r_date_match[1];
            sig = r_date_match[2];
            date_name.push_back((std::string) "");
            c14_age.push_back(std::strtod(age.c_str(), nullptr));
            c14_sig.push_back(std::strtod(sig.c_str(), nullptr));
        } else if (np_model && regex_search(line, r_date_match, named_r_f14c_regex)) {
            name = r_date_match[1];
            age = r_date_match[2];
            sig = r_date_match[3];
            date_name.push_back(name);
            f14c_age.push_back(std::strtod(age.c_str(), nullptr));
            f14c_sig.push_back(std::strtod(sig.c_str(), nullptr));
        } else if (np_model && regex_search(line, r_date_match, unnamed_r_f14c_regex)){
            age = r_date_match[1];
            sig = r_date_match[2];
            date_name.push_back((std::string) "");
            f14c_age.push_back(std::strtod(age.c_str(), nullptr));
            f14c_sig.push_back(std::strtod(sig.c_str(), nullptr));
        }
    }
    if (!c14_age.empty() && !f14c_age.empty())  throw InconsistentDateFormatsException();
    return np_model && (!c14_age.empty() || !f14c_age.empty());
}

/* In case multiple different models are specified in OxCal, we need to determine the offset of the output index
 * for the NP model, so that when we produce output the indices match those already given by OxCal in the output
 * file.
 */
int read_output_offset(const std::string& model_name) {
    std::string line, model_index;
    std::regex np_output_regex(R"(ocd\[([0-9]+)\].name\s*=\s*["'])" + model_name + R"(["'];)");
    std::smatch np_model_match;

    std::string output_file_path = project_directory + project_name + ".js";
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


void read_oxcal_version() {
    std::string line, model_index;
    std::regex oxcal_version_regex(R"(ocd\[0\].ref\s*=\s*["'](.*)["'];)");
    std::smatch oxcal_version_match;

    std::string output_file_path = project_directory + project_name + ".js";
    std::fstream file(output_file_path, std::ios::in);
    if(!file.is_open()) throw UnableToReadOutputFileException(output_file_path);

    oxcal_version = "OxCal, Bronk Ramsey (2021)"; // Used if we cannot find the version in the output
    while (getline(file, line)) {
        if (regex_search(line, oxcal_version_match, oxcal_version_regex)) {
            oxcal_version = oxcal_version_match[1];
        }
    }
}


/* Reads in the options from the oxcal data file. If any of the options are not found in the file
 * then the value will not be altered from the original value. Currently, it will populate the
 * following option variables provided as arguments:
 * - iterations: The number of iterations for the DPMM
 * - resolution: The resolution used for outputting the predictive and posterior density
 * - ranges: A vector of 3 values denoting whether to log the 68.3%, 95.4% and 99.7% ranges
 * - quantile ranges: False to use HPD ranges, True to simply work out a single quantile range
 * - calibration_curve_name: Which calibration curve to use
 */
void read_default_options_from_data_file(
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        std::string &calibration_curve_name
        ) {
    std::string filepath = calling_directory + "OxCal.dat";
    std::fstream file(filepath, std::ios::in);
    if(!file.is_open()) throw UnableToReadDefaultOptionsFileException(filepath);

    std::regex option_regex(R"(-(\w)(.*))");
    std::smatch option_match;
    std::string line, option, value;

    while (getline(file, line)) {
        if (regex_search(line, option_match, option_regex)) {
            option = option_match[1];
            value = option_match[2];
            if (option == "i") {
                resolution = std::stod(value);
            } else if (option == "s") {
                ranges[value[0] - 49] = value[1] == '1'; // Subtract by 49 to convert ascii value for 1, 2, 3 to integer
            } else if (option == "h") {
                quantile_ranges = value == "1";
            } else if (option == "c") {
                calibration_curve_name = value;
            }
        }
    }
}

/* Reads in the options from the oxcal data file. If any of the options are not found in the file
 * then the value will not be altered from the original value. Currently, it will populate the
 * following option variables provided as arguments:
 * - iterations: The number of iterations for the DPMM
 * - resolution: The resolution used for outputting the predictive and posterior density
 * - ranges: A vector of 3 values denoting whether to log the 68.3%, 95.4% and 99.7% ranges
 * - quantile ranges: False to use HPD ranges, True to simply work out a single quantile range
 * - use_f14c: Whether to perform the calculations in f14c space
 * - calibration_curve_name: Which calibration curve to use
 */
void read_options_from_oxcal_file(
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        bool &use_f14c,
        std::string &calibration_curve_name,
        int &seed) {

    std::string line, option, value, end_of_section = "};";
    std::regex options_regex(R"(Options\(\s*\))");
    std::regex option_regex(R"((\w+)=['"]*([^'"]+)['"]*;)");
    bool options_block = false;
    std::smatch option_match, m;
    std::regex allowed_calibration_curve_regex(R"((intcal|shcal)[0-9]+\.14c)");

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
                if (!regex_search(value, m, allowed_calibration_curve_regex)) {
                    std::string log_string = "The calibration curve " + value;
                    log_string += " is not recognised.\nDefault " + calibration_curve_name + " is being used.\n";
                    update_log_file(log_string);
                } else {
                    calibration_curve_name = value;
                }
            } else if (option == "RandomNumberSeed") {
                seed = std::stoi(value);
            }
        }
    }
}
