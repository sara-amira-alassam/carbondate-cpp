#include <iostream>
#include "src/WalkerDPMM.h"
#include "src/csv_helpers.h"

int main() {
    std::vector<double> c14_age = get_csv_data_from_column("../data/kerr.csv", 0);
    std::vector<double> c14_sig = get_csv_data_from_column("../data/kerr.csv", 1);
    std::vector<double> cc_cal_age = get_csv_data_from_column("../data/intcal20.14c", 0);
    std::vector<double> cc_c14_age = get_csv_data_from_column("../data/intcal20.14c", 1);
    std::vector<double> cc_c14_sig = get_csv_data_from_column("../data/intcal20.14c", 2);
    WalkerDPMM dpmm;

    dpmm.initialise(c14_age, c14_sig, cc_cal_age, cc_c14_age, cc_c14_sig);
    dpmm.calibrate(1e5, 10);
    int ident = 10;
    write_column_to_csv(
            "../output/kerr_calendar_age_" + std::to_string(ident) + ".csv",
            "calendar age",
            dpmm.get_calendar_age(ident));

    DensityData predictive_density = dpmm.get_predictive_density(500, 101, 0.025);
    std::vector<std::vector<double>> density_data = {
        predictive_density.cal_age,
        predictive_density.ci_lower,
        predictive_density.mean,
        predictive_density.ci_upper
    };
    std::vector<std::string> density_headers = {"calendar_age", "ci_lower", "mean", "ci_upper"};
    write_columns_to_csv("../output/kerr_predictive_density.csv", density_headers, density_data);

    DensityOutput density_output = dpmm.get_posterior_density(1, ident, 5);
    density_output.print();

    std::vector<std::string> post_density_headers = {"calendar_age", "probability"};
    write_columns_to_csv(
            "../output/kerr_posterior_density_" + std::to_string(ident) + ".csv",
            post_density_headers,
            density_output.as_columns());

    return 0;
}
