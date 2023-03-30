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
        bool quantile_ranges,
        const std::vector<bool>& log_ranges,
        const std::vector<double>& posterior_calendar_ages_AD)
        : _log_ranges(log_ranges), DensityOutput(ident + offset + 1, resolution) {

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

    if (quantile_ranges) {
        double edge_width;
        for (double probability : _range_probabilities) {
            edge_width = (1. - probability) / 2.;
            std::vector<std::vector<double>> range(1, std::vector<double> {0, 0, probability});
            std::vector<double> calendar_ages(posterior_calendar_ages_AD);
            edge_quantiles(calendar_ages, edge_width, range[0][0], range[0][1]);
            _ranges.push_back(range);
        }
    } else {
        for (double probability : _range_probabilities) {
            _ranges.push_back(get_ranges_by_bisection(probability));
        }
    }
}

std::string PosteriorDensityOutput::range_lines(int& comment_index) {
    std::string range_lines;

    for (int i = 0; i < _ranges.size(); i++) {
        std::string range_string = "range[" + std::to_string(i + 1) + "]";
        if (i == 0) {
            range_lines = _output_prefix + ".range=[];\n";
        }
        range_lines += _output_prefix + "." + range_string + "=[];\n";
        for (int j = 0; j < _ranges[i].size(); j++) {
            range_lines += output_line(range_string + "[" + std::to_string(j) + "]", _ranges[i][j]);
        }
        if (_log_ranges[i]) {
            range_lines += comment_line(
                    "  " + to_percent_string(_range_probabilities[i]) + " probability", comment_index);
            for (auto & range : _ranges[i]) {
                std::string comment = "    " + std::to_string(int (round(range[0]))) + "AD";
                comment += " (" + to_percent_string(range[2] / 100.) + ") ";
                comment += std::to_string(int (round(range[1]))) + "AD";
                range_lines += comment_line(comment, comment_index);
            }
        }
    }
    return range_lines;
}
