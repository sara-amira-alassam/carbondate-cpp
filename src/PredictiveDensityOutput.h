//
// Created by Sara Al-Assam on 20/03/2023.
//

#ifndef CARBONDATE_PREDICTIVEDENSITYOUTPUT_H
#define CARBONDATE_PREDICTIVEDENSITYOUTPUT_H
#include "DensityOutput.h"

class PredictiveDensityOutput : public DensityOutput {
    std::vector<double> _ci_lower;
    std::vector<double> _ci_upper;
    int _n_obs;

private:
    std::vector<std::string> get_output_lines() override;

public:
    PredictiveDensityOutput(int n_obs, int offset, double resolution);
    void set_confidence_intervals(const std::vector<double>& ci_lower, const std::vector<double>& ci_upper);
};


#endif //CARBONDATE_PREDICTIVEDENSITYOUTPUT_H
