#include <iostream>
#include "src/WalkerDPMM.h"
#include "src/csv_helpers.h"
#include "src/read_data.h"

int main(int argc, char* argv[]) {
    if (argc < 2)
        return 1;
    std::string file_prefix = argv[1];

    std::vector<double> c14_age, c14_sig;
    std::vector<std::string> c14_name;
    std::vector<double> cc_cal_age, cc_c14_age, cc_c14_sig;
    WalkerDPMM dpmm;

    read_oxcal_data(file_prefix, c14_name, c14_age, c14_sig);
    read_calibration_curve("../data/intcal20.14c", cc_cal_age, cc_c14_age, cc_c14_sig);

    dpmm.initialise(c14_age, c14_sig, c14_name, cc_cal_age, cc_c14_age, cc_c14_sig);

    for (int i = 0; i <= 10; i++){
        DensityOutput density_output = dpmm.get_single_calendar_age_likelihood(i);
        density_output.print(5);
        density_output.write_to_file(5, file_prefix);
    }

    dpmm.calibrate(1e5, 10);

    for (int i = 0; i <= 10; i++){
        DensityOutput density_output = dpmm.get_posterior_calendar_age_density(i);
        density_output.print(5);
        density_output.write_to_file(5, file_prefix);
    }

    /*
    DensityData predictive_density = dpmm.get_predictive_density(5000, 101, 0.025);
    std::vector<std::vector<double>> density_data = {
        predictive_density.cal_age,
        predictive_density.ci_lower,
        predictive_density.mean,
        predictive_density.ci_upper
    };
    std::vector<std::string> density_headers = {"calendar_age", "ci_lower", "mean", "ci_upper"};
    write_columns_to_csv("../output/kerr_predictive_density.csv", density_headers, density_data);*/

    return 0;
}
