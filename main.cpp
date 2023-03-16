#include "src/WalkerDPMM.h"
#include "src/OxCalOutput.h"
#include "src/read_data.h"

int main(int argc, char* argv[]) {
    if (argc < 2)
        return 1;
    std::string file_prefix = argv[1];

    const int output_resolution = 5;
    std::vector<double> c14_age, c14_sig;
    std::vector<std::string> c14_name;
    std::vector<double> cc_cal_age, cc_c14_age, cc_c14_sig;
    WalkerDPMM dpmm;

    if (!read_oxcal_data(file_prefix, c14_name, c14_age, c14_sig)) {
        // If there is no data within the NP model in this OxCal file then simply exit
        return 0;
    }
    read_calibration_curve("../data/intcal20.14c", cc_cal_age, cc_c14_age, cc_c14_sig);
    OxCalOutput oxcal_output(11, output_resolution, file_prefix);

    dpmm.initialise(c14_age, c14_sig, c14_name, cc_cal_age, cc_c14_age, cc_c14_sig);
    dpmm.calibrate(5e5, 10);

    for (int i = 0; i < c14_age.size(); i++){
        oxcal_output.set_posterior(i, dpmm.get_posterior_calendar_age_density(i));
    }

    oxcal_output.print();

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
