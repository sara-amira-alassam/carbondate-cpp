//
// Created by Sara Al-Assam on 20/03/2023.
//

#ifndef CARBONDATE_POSTERIORDENSITYOUTPUT_H
#define CARBONDATE_POSTERIORDENSITYOUTPUT_H

#include "DensityOutput.h"

class PosteriorDensityOutput : public DensityOutput {
    const std::vector<double> _range_probabilities{0.68268949, 0.95449974, 0.99730020};
    // Top level vector index represents each probability
    // Next level vector index is the section of the range
    // Final level vector index gives start point, end point, and probability of that section
    std::vector<std::vector<std::vector<double>>> _ranges;
    std::vector<bool> _log_ranges;
    bool _intercept_method;

    std::string range_lines(int &comment_index) override;

public:
    PosteriorDensityOutput(int ident, int offset, double resolution,
                           bool quantile_ranges, bool intercept_ranges,
                           const std::vector<bool> &log_ranges,
                           const std::vector<double> &posterior_calendar_ages_AD);
};


#endif //CARBONDATE_POSTERIORDENSITYOUTPUT_H
