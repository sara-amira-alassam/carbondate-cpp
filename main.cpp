#include "src/WalkerDPMM.h"
#include "src/PolyaUrnDPMM.h"
#include "src/read_data.h"
#include "src/PredictiveDensityOutput.h"
#include "src/PosteriorDensityOutput.h"


int main(int argc, char* argv[]) {
    if (argc < 2)
        return 1;
    std::string file_prefix = argv[1];

    const int n_posterior_samples = 5000;
    const double quantile_edge_width = 0.1586553; // 1-sigma interval
    int output_offset;
    std::vector<double> c14_age, c14_sig, f14c_age, f14c_sig;
    std::string model_name, calibration_curve = "intcal20.14c";
    std::vector<double> cc_cal_age, cc_c14_age, cc_c14_sig;
    PolyaUrnDPMM dpmm;

    // The following relate to options that may be overwritten in the call below to read_options()
    int num_iterations = 1e5;
    double output_resolution = 5;
    std::vector<bool> log_ranges {true, true, false}; // log 1, 2, 3 s.d. ranges respectively?
    bool quantile_ranges = false, use_f14c = true;

    if (!read_oxcal_data(file_prefix, c14_age, c14_sig, f14c_age, f14c_sig, model_name)) {
        // If there is no data within the NP model in this OxCal file then simply exit
        return 0;
    }
    output_offset = read_output_offset(file_prefix, model_name);
    read_options(
        file_prefix,
        num_iterations,
        output_resolution,
        log_ranges,
        quantile_ranges,
        use_f14c,
        calibration_curve);
    read_calibration_curve(calibration_curve, cc_cal_age, cc_c14_age, cc_c14_sig);

    if (use_f14c) {
        if (f14c_age.empty()) convert_to_f14c_age(c14_age, c14_sig, f14c_age, f14c_sig);
        dpmm.initialise(f14c_age, f14c_sig, true, cc_cal_age, cc_c14_age, cc_c14_sig, 0);
    } else {
        if (c14_age.empty()) convert_to_c14_age(f14c_age, f14c_sig, c14_age, c14_sig);
        dpmm.initialise(c14_age, c14_sig, false, cc_cal_age, cc_c14_age, cc_c14_sig, 0);
    }
    dpmm.calibrate(num_iterations, 10);

    DensityData predictive_density_data = dpmm.get_predictive_density(
            n_posterior_samples, output_resolution, quantile_edge_width);
    PredictiveDensityOutput predictive_density(
            dpmm.get_nobs(),
            output_offset,
            output_resolution,
            model_name,
            predictive_density_data.cal_age_AD,
            predictive_density_data.mean,
            predictive_density_data.ci_lower,
            predictive_density_data.ci_upper);
    predictive_density.print(file_prefix);

    for (int i = 0; i < dpmm.get_nobs(); i++){
        PosteriorDensityOutput posterior_density(
                i,
                output_offset,
                output_resolution,
                quantile_ranges,
                log_ranges,
                dpmm.get_posterior_calendar_ages(i));
        posterior_density.print(file_prefix);
    }

    return 0;
}
