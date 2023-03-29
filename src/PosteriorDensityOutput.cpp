//
// Created by Sara Al-Assam on 20/03/2023.
//
#include <cmath>
#include "helpers.h"
#include "PosteriorDensityOutput.h"

// Creates an object suitable for printing out the posterior calendar age density for each
// determination, taking the following arguments:
// * ident:  The index for the determination corresponding to this posterior
// * offset: The offset to use for indexing the total model output
//           (e.g. if other models have been specified in the input file in addition to this one)
// * resolution: The resolution to use for detmining the probability curve
// * posterior_calendar_ages: A list of sampled calendar ages for this determination, in calAD
PosteriorDensityOutput::PosteriorDensityOutput(
        int ident,
        int offset,
        double resolution,
        const std::vector<bool>& ranges,
        const std::vector<double>& posterior_calendar_ages_AD)
        : DensityOutput(ident + offset + 1, resolution, ranges) {

    unsigned n = posterior_calendar_ages_AD.size();
    double min_calendar_age = std::numeric_limits<double>::infinity();
    double max_calendar_age = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        if (posterior_calendar_ages_AD[i] < min_calendar_age) {
            min_calendar_age = posterior_calendar_ages_AD[i];
        } else if (posterior_calendar_ages_AD[i] > max_calendar_age) {
            max_calendar_age = posterior_calendar_ages_AD[i];
        }
    }
    _start_calAD = floor(min_calendar_age + resolution / 2.0);

    int num_breaks = (int) ceil((max_calendar_age - min_calendar_age) / resolution) + 1;
    double break_start = _start_calAD - resolution / 2.0;
    std::vector<double> _probability(num_breaks, 0);
    for (int i = 0; i < n; i++) {
        _probability[(int) ((posterior_calendar_ages_AD[i] - break_start) / resolution)]++;
    }
    set_probability(_probability);
    _mean_calAD = mean(posterior_calendar_ages_AD);
    _median_calAD = median(posterior_calendar_ages_AD);
    _sigma_calAD = sigma(posterior_calendar_ages_AD, _mean_calAD);
}
