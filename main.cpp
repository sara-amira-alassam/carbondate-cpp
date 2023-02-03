#include <iostream>
#include "src/WalkerDPMM.h"

std::vector<double> get_csv_data_from_column(const std::string& filename,int column_index);

int main() {
    std::vector<double> c14_age = get_csv_data_from_column("../data/kerr.csv", 2);
    std::vector<double> c14_sig = get_csv_data_from_column("../data/kerr.csv", 3);
    std::vector<double> cc_cal_age = get_csv_data_from_column("../data/intcal20.14c", 0);
    std::vector<double> cc_c14_age = get_csv_data_from_column("../data/intcal20.14c", 1);
    std::vector<double> cc_c14_sig = get_csv_data_from_column("../data/intcal20.14c", 2);
    WalkerDPMM dpmm;

    dpmm.initialise(c14_age, c14_sig, cc_cal_age, cc_c14_age, cc_c14_sig);
    dpmm.calibrate(10000, 10);
    dpmm.get_predictive_density(5000, 1001, 0.025);

    printf("dpmm %f\n", dpmm.get_alpha()[0]);

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
