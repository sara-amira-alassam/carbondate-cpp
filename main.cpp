#include <iostream>
#include "src/WalkerDPMM.h"

int main() {
    std::vector<double> c14_age{1421, 1551, 1475, 1289, 1270, 1552, 1191, 1475, 1585, 1470};
    std::vector<double> c14_sig{32, 33, 21, 32, 21, 49, 30, 60, 70, 120};
    std::vector<double> cc_cal_age{1000, 1200, 1400, 1600, 1800};
    std::vector<double> cc_c14_age{1126, 1301, 1537, 1714, 1887};
    std::vector<double> cc_c14_sig{13, 10, 13, 12, 14};

    WalkerDPMM dpmm;

    dpmm.initialise(c14_age, c14_sig, cc_cal_age, cc_c14_age, cc_c14_sig);

    printf("dpmm %f\n", dpmm.get_alpha()[0]);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
