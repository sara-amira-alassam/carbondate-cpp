#ifndef CARBONDATE_READ_DATA_H
#define CARBONDATE_READ_DATA_H

#include "carbondate.h"

class IncorrectArgumentsException : public CarbondateException {
private:
    std::string _error_message = "Incorrect arguments - correct usage: carbondate [project_name]";
};

class UnableToReadDefaultOptionsFileException : public CarbondateException {
public:
    explicit UnableToReadDefaultOptionsFileException(const std::string& file_path) {
        _error_message = "Unable to read file with default options at " + file_path;
    }
};

class UnableToReadCalibrationCurveException : public CarbondateException {
public:
    explicit UnableToReadCalibrationCurveException(const std::string& calibration_curve_name) {
        _error_message = "Unable to read calibration curve " + calibration_curve_name;
    }
};

class UnableToReadOxcalFileException : public CarbondateException {
public:
    explicit UnableToReadOxcalFileException(const std::string& file_path) {
        _error_message = "Unable to read oxcal file " + file_path;
    }
};

class UnableToReadOutputFileException : public CarbondateException {
public:
    explicit UnableToReadOutputFileException(const std::string& file_path) {
        _error_message = "Unable to read output file " + file_path;
    }
};

class UnableToDetermineOutputOffsetException : public CarbondateException {
public:
    explicit UnableToDetermineOutputOffsetException(const std::string& file_path, const std::string& model_name) {
        _error_message = "Unable to find NP model `" + model_name + "' in output file " + file_path;
    }
};

class InconsistentDateFormatsException : public CarbondateException {
public:
    explicit InconsistentDateFormatsException() {
        _error_message = "All dates must be the same format.";
    }
};

void read_arguments(int argc, char* argv[]);

void read_calibration_curve(
    const std::string& calibration_curve_name,
    std::vector<double>& cc_cal_age,
    std::vector<double>& cc_c14_age,
    std::vector<double>& cc_c14_sig);

bool read_oxcal_data(
    std::vector<std::string>& date_name,
    std::vector<double>& c14_age,
    std::vector<double>& c14_sig,
    std::vector<double>& f14c_age,
    std::vector<double>& f14c_sig,
    std::string& model_name);

int read_output_offset(const std::string& model_name);

void read_default_options_from_data_file(
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        std::string &calibration_curve_name
);

void read_options_from_oxcal_file(
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_ranges,
        bool &use_f14c,
        std::string &calibration_curve_name);

#endif //CARBONDATE_READ_DATA_H
