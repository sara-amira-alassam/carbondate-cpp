#include "src/WalkerDPMM.h"
#include "src/OxCalOutput.h"
#include "src/read_data.h"
#include "src/csv_helpers.h"

int main(int argc, char* argv[]) {
    if (argc < 2)
        return 1;
    std::string file_prefix = argv[1];

    const int output_resolution = 5, output_offset = 0, n_posterior_samples = 5000;
    const double quantile_edge_width = 0.1586553; // 1-sigma interval
    std::vector<double> c14_age, c14_sig;
    std::string model_name;
    std::vector<double> cc_cal_age, cc_c14_age, cc_c14_sig;
    WalkerDPMM dpmm;

    if (!read_oxcal_data(file_prefix, c14_age, c14_sig, model_name)) {
        // If there is no data within the NP model in this OxCal file then simply exit
        return 0;
    }
    read_calibration_curve("../data/intcal20.14c", cc_cal_age, cc_c14_age, cc_c14_sig);
    OxCalOutput oxcal_output(file_prefix);

    dpmm.initialise(c14_age, c14_sig, cc_cal_age, cc_c14_age, cc_c14_sig);
    dpmm.calibrate(1e5, 10);

    DensityData predictive_density_data = dpmm.get_predictive_density(
            n_posterior_samples, output_resolution, quantile_edge_width);
    std::vector<std::string> headers = {"date", "mean", "ci_lower", "ci_upper"};
    std::vector<std::vector<double>> data = {
            predictive_density_data.cal_age_AD,
            predictive_density_data.mean,
            predictive_density_data.ci_lower,
            predictive_density_data.ci_upper
    };
    write_columns_to_csv("../output/np_kerr.csv", headers, data);

    PredictiveDensityOutput predictive_density(
            (int) c14_age.size(),
            output_offset,
            output_resolution,
            model_name,
            predictive_density_data.cal_age_AD,
            predictive_density_data.mean,
            predictive_density_data.ci_lower,
            predictive_density_data.ci_upper);
    //predictive_density.print(file_prefix);

    for (int i = 0; i < c14_age.size(); i++){
        PosteriorDensityOutput posterior_density(
                i, output_offset, output_resolution, dpmm.get_posterior_calendar_ages(i));
        //posterior_density.print(file_prefix);
    }

    return 0;
}
