/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <cmath>
#include "PosteriorDensityOutput.h"

/*
 * Creates an object suitable for printing out the posterior calendar age density for each
 * determination, taking the following arguments:
 * - ident:  The index for the determination corresponding to this posterior
 * - date_name: The name given to the data point, or an empty string if no name is given
 * - rc_age: The radiocarbon determination
 * - rc_sig: The error for the radiocarbon determination
 * - f14c_age: True if the radiocarbon determination is F14C age, false otherwise
 * - offset: The offset to use for indexing the total model output
 *          (e.g. if other models have been specified in the input file in addition to this one)
 * - resolution:
 *    The resolution to use for determining the probability curve. Always rounded up to the
 *    nearest whole number.
 * - quantile_ranges: True if quantile ranges should be used, false otherwise
 * - log_ranges: A vector of length 3, denoting whether the ranges should be logged for each of the range probabilities
 * - posterior_calendar_ages: A list of sampled calendar ages for this determination, in calAD
 */
PosteriorDensityOutput::PosteriorDensityOutput(
        int ident,
        const std::string& date_name,
        double rc_age,
        double rc_sig,
        bool f14c_age,
        int offset,
        double resolution,
        bool quantile_ranges,
        const std::vector<bool> &log_ranges,
        const std::vector<double> &posterior_calendar_ages_AD)
        : _log_ranges(log_ranges), DensityOutput(ident + offset + 1, resolution) {

    if (date_name.length() > 0)
        _label = date_name;
    else if (f14c_age)
        _label = "R_F14C(" + to_string(rc_age, 6) + "," + to_string(rc_sig, 8) + ")";
    else
        _label = "R_Date(" + to_string(rc_age, 4) + "," + to_string(rc_sig, 4) + ")";

    unsigned n = posterior_calendar_ages_AD.size();
    double min_calendar_age = posterior_calendar_ages_AD[0];
    double max_calendar_age = posterior_calendar_ages_AD[0];
    for (int i = 1; i < n; i++) {
        if (posterior_calendar_ages_AD[i] < min_calendar_age) {
            min_calendar_age = posterior_calendar_ages_AD[i];
        } else if (posterior_calendar_ages_AD[i] > max_calendar_age) {
            max_calendar_age = posterior_calendar_ages_AD[i];
        }
    }
    _start_calAD = floor(min_calendar_age - resolution / 2.0);

    double break_start = _start_calAD - resolution / 2.0;
    int num_breaks = (int) ceil((max_calendar_age - break_start) / resolution);
    std::vector<double> probability(num_breaks, 0);
    for (int i = 0; i < n; i++) probability[(int) ((posterior_calendar_ages_AD[i] - break_start) / resolution)]++;
    _set_probability(probability);

    _mean_calAD = mean(posterior_calendar_ages_AD);
    _median_calAD = median(posterior_calendar_ages_AD);
    _sigma_calAD = sigma(posterior_calendar_ages_AD, _mean_calAD);

    if (quantile_ranges) {
        double edge_width;
        for (double range_probability : _range_probabilities) {
            edge_width = (1. - range_probability) / 2.;
            std::vector<std::vector<double>> range(1, std::vector<double> {0, 0, range_probability * 100.});
            std::vector<double> calendar_ages(posterior_calendar_ages_AD);
            edge_quantiles(calendar_ages, edge_width, range[0][0], range[0][1]);
            _ranges.push_back(range);
        }
    } else {
        for (double range_probability : _range_probabilities) {
            _ranges.push_back(get_ranges_by_intercepts(range_probability));
        }
    }
}

// Returns the area under the probability curve if we ignore all values below the cut-off
// Also populates the vector ranges, where each entry contains
// [start_calAD, end_calAD, probability within this range]
double PosteriorDensityOutput::find_probability_and_ranges_for_cut_off(
        double cut_off, std::vector<std::vector<double>> &ranges) {
    ranges.clear();
    double y1, y2, dx, x_intercept_1, x_intercept_2, res = _resolution;
    double range_probability = 0, total_probability = 0;
    const double min_prob = 0.001; // Don't bother to store ranges with probability less than this
    for (int i = 0; i < _probability.size() - 1; i++) {
        y1 = _probability[i];
        y2 = _probability[i + 1];
        if (i == 0 && y1 > cut_off) {
            x_intercept_1 = _start_calAD + i * res;
            range_probability = (y1 + y2) * res / 2.;
        } else if (y1 <= cut_off and y2 > cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_1 = _start_calAD + i * res + dx;
            range_probability = (cut_off + y2) * (res - dx) / 2.;
        } else if (y1 > cut_off and y2 <= cut_off) {
            dx = res * (cut_off - y1) / (y2 - y1);
            x_intercept_2 = _start_calAD + i * res + dx;
            range_probability += (y1 + cut_off) * dx / 2.;
            range_probability *= _prob_norm;
            if (range_probability > min_prob) {
                ranges.push_back(std::vector<double> {x_intercept_1, x_intercept_2, range_probability});
                total_probability += range_probability;
            }
            range_probability = 0;
        } else if (y1 > cut_off and y2 > cut_off) {
            range_probability += (y1 + y2) * res / 2.;
        }
        if (y2 > cut_off and i == _probability.size() - 2) {
            x_intercept_2 = _start_calAD + (i + 1) * res;
            range_probability *= _prob_norm;
            if (range_probability > min_prob) {
                ranges.push_back(std::vector<double> {x_intercept_1, x_intercept_2, range_probability});
                total_probability += range_probability;
            }
        }
    }
    for (std::vector<double> &range : ranges) range[2] *= 100.;
    return total_probability;
}

// Merges any ranges which are separated by a distance less than the current resolution
std::vector<std::vector<double>> PosteriorDensityOutput::simplify_ranges(std::vector<std::vector<double>> &ranges) {
    std::vector<std::vector<double>> new_ranges;

    new_ranges.push_back(ranges[0]);
    for (int i = 1; i < ranges.size(); i++) {
        if (ranges[i][0] - ranges[i - 1][1] >= _resolution) {
            new_ranges.push_back(ranges[i]);
        } else {
            new_ranges[new_ranges.size() - 1][1] = ranges[i][1];
            new_ranges[new_ranges.size() - 1][2] += ranges[i][2];
        }
    }
    return new_ranges;
}

// Finds the calendar age ranges between which the probability matches the provided probability
// where the points with the highest probability are chosen first. Returns the result as a
// vector of vectors, where each entry contains
// [start_calAD, end_calAD, probability within this range]
// arg `resolution` gives the resolution of the probability curve to use
std::vector<std::vector<double>> PosteriorDensityOutput::get_ranges_by_intercepts(double probability) {
    std::vector<std::vector<double>> ranges;

    // Use bisection method to find the closest probability cut-off to give the desired probability
    // a and b are the upper and lower points of the section - we know the smoothed probability has a max of 1
    double a = 0., b = 1.;
    double p; // p is the midpoint between a and b
    double current_probability;
    const int max_iter = 1000;
    for (int i = 0; i < max_iter; i++) {
        p = (a + b) / 2.;
        current_probability = find_probability_and_ranges_for_cut_off(p, ranges);
        if (fabs(current_probability - probability) < 1e-4) {
            break;
        }
        if (current_probability < probability) {
            b = p;
        } else {
            a = p;
        }
    }
    ranges = simplify_ranges(ranges);
    return ranges;
}

std::string PosteriorDensityOutput::_range_lines(int &comment_index) {
    std::string range_lines, comment, log_lines = "Posterior " + _label + "\n";
    std::vector<double> text_ranges;

    for (int i = 0; i < _ranges.size(); i++) {
        std::string range_string = "range[" + std::to_string(i + 1) + "]";
        if (i == 0) {
            range_lines = _output_prefix + ".range=[];\n";
        }
        range_lines += _output_prefix + "." + range_string + "=[];\n";
        for (int j = 0; j < _ranges[i].size(); j++) {
            range_lines += _output_line(range_string + "[" + std::to_string(j) + "]", _ranges[i][j]);
        }
        if (_log_ranges[i]) {
            comment = "  " + to_percent_string(_range_probabilities[i]) + " probability";
            range_lines += _comment_line(comment, comment_index);
            log_lines += comment + "\n";
            for (std::vector<double> & range : _ranges[i]) {
                comment = "    " + std::to_string(int (round(range[0]))) + "AD";
                comment += " (" + to_percent_string(range[2] / 100.) + ") ";
                comment += std::to_string(int (round(range[1]))) + "AD";
                range_lines += _comment_line(comment, comment_index);
                log_lines += comment + "\n";
            }
            text_ranges.push_back(_ranges[i][0][0]);
            text_ranges.push_back(_ranges[i][_ranges[i].size() - 1][1]);
        }
    }
    _log_lines.push_back(log_lines.substr(0, log_lines.length() - 2)); // Remove the last carriage return
    if (!text_ranges.empty()) _add_text_ranges(text_ranges);
    return range_lines;
}


void PosteriorDensityOutput::_add_text_ranges(const std::vector<double> &values) {
    std::string text_line = "@" + _label;

    for (double value : values) text_line += "\t" + to_string(value, 4);

    _text_lines.push_back(text_line);
}
