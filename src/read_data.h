#ifndef CARBONDATE_READ_DATA_H
#define CARBONDATE_READ_DATA_H

#include "carbondate.h"

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

void read_calibration_curve(
    const std::string& calibration_curve,
    std::vector<double>& cc_cal_age,
    std::vector<double>& cc_c14_age,
    std::vector<double>& cc_c14_sig);

bool read_oxcal_data(
    const std::string& file_prefix,
    std::vector<double>& c14_age,
    std::vector<double>& c14_sig,
    std::vector<double>& f14c_age,
    std::vector<double>& f14c_sig,
    std::string& model_name);

int read_output_offset(const std::string& file_prefix, const std::string& model_name);

void read_options(
        const std::string &file_prefix,
        int &iterations,
        double &resolution,
        std::vector<bool> &ranges,
        bool &quantile_range,
        bool &use_f14c,
        std::string &calibration_curve);

#endif //CARBONDATE_READ_DATA_H
