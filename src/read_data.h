#ifndef CARBONDATE_READ_DATA_H
#define CARBONDATE_READ_DATA_H

#include <string>
#include <vector>

#endif //CARBONDATE_READ_DATA_H

void read_calibration_curve(
    const std::string& calibration_curve_path,
    std::vector<double>& cc_cal_age,
    std::vector<double>& cc_c14_age,
    std::vector<double>& cc_c14_sig);

void read_oxcal_data(
        const std::string& file_prefix,
        std::vector<std::string>& c14_name,
        std::vector<double>& c14_age,
        std::vector<double>& c14_sig);
