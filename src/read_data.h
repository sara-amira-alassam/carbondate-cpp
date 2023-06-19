#ifndef CARBONDATE_READ_DATA_H
#define CARBONDATE_READ_DATA_H

#include <string>
#include <vector>

#endif //CARBONDATE_READ_DATA_H

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
        bool &intercept_ranges,
        std::string &calibration_curve);
