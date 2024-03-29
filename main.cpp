/* carbondate Copyright (C) 2004 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <iostream>
#include "inc/carbondate.h"

int main(int argc, char* argv[]) {

    const int n_posterior_samples = 5000;
    const double quantile_edge_width = 0.1586553; // 1-sigma interval
    int output_offset;  // This takes into account any results calculated by Oxcal before the NP model data
    std::vector<double> c14_age, c14_sig, f14c_age, f14c_sig;  // Input data - radiocarbon ages and errors
    std::string model_name, calibration_curve = "intcal20.14c";
    std::vector<double> cc_cal_age, cc_c14_age, cc_c14_sig;
    std::vector<std::string> date_name;
    std::vector<std::string> js_output_lines, log_lines, text_lines;
    PolyaUrnDPMM dpmm;

    // The following relate to options that are read in from the OxCal.dat file, and also may be overwritten in the
    // call below to read_options_from_oxcal_file()
    int num_iterations = 1e5;
    double output_resolution;
    std::vector<bool> log_ranges(3); // log 1, 2, 3 s.d. ranges respectively?
    bool quantile_ranges, use_f14c = true; // use_f14c is the only hard-coded option

    // This option is for the random number seed, it will not be documented so only used by developers.
    // The default is zero (so the seed is chosen based on the time i.e. different for every run), but it can
    // be set to a non-zero integer for reproducible results.
    int seed = 0;

    // How much to thin the samples
    int n_thin = 10;

    std::cout << carbondate_full_reference() << std::endl;

    try {
        read_arguments(argc, argv);
        create_work_file();
        initialize_log_file();

        if (!read_oxcal_data(date_name, c14_age, c14_sig, f14c_age, f14c_sig, model_name)) {
            // If there is no data within the NP model in this OxCal file then simply exit
            return 0;
        }
        output_offset = read_output_offset(model_name);
        read_default_options_from_data_file(output_resolution, log_ranges, quantile_ranges, calibration_curve);
        read_options_from_oxcal_file(
                num_iterations, output_resolution, log_ranges, quantile_ranges, use_f14c, calibration_curve, seed);
        read_calibration_curve(calibration_curve, cc_cal_age, cc_c14_age, cc_c14_sig);
        read_oxcal_version();

        if (use_f14c) {
            if (f14c_age.empty()) convert_to_f14c_age(c14_age, c14_sig, f14c_age, f14c_sig);
            dpmm.initialise(f14c_age, f14c_sig, true, cc_cal_age, cc_c14_age, cc_c14_sig, seed);
        } else {
            if (c14_age.empty()) convert_to_c14_age(f14c_age, f14c_sig, c14_age, c14_sig);
            dpmm.initialise(c14_age, c14_sig, false, cc_cal_age, cc_c14_age, cc_c14_sig, seed);
        }
        num_iterations = dpmm.calibrate(num_iterations, n_thin);

        update_work_file_postprocessing(num_iterations);

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
        predictive_density.append_output(js_output_lines, log_lines, text_lines);

        for (int i = 0; i < dpmm.get_nobs(); i++){
            PosteriorDensityOutput posterior_density(
                    i,
                    date_name[i],
                    c14_age.empty() ? f14c_age[i] : c14_age[i],
                    c14_age.empty() ? f14c_sig[i] : c14_sig[i],
                    c14_age.empty(),
                    output_offset,
                    output_resolution,
                    quantile_ranges,
                    log_ranges,
                    dpmm.get_posterior_calendar_ages(i));
            posterior_density.append_output(js_output_lines, log_lines, text_lines);
        }

        update_js_output_file(js_output_lines);
        update_log_file(log_lines);
        update_text_file(text_lines);

    } catch (const CarbondateException& ex) {
        std::cout << "Exception caught: " << ex.what() << std::endl;
        remove_work_file();
        return 1;
    }

    remove_work_file();
    return 0;
}
