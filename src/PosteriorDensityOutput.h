/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
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
    std::string _label;

    std::string _range_lines(int &comment_index) override;

public:
    PosteriorDensityOutput(
            int ident, const std::string& date_name, double rc_age, double rc_sig, bool f14c_age, int offset,
            double resolution, bool quantile_ranges, const std::vector<bool> &log_ranges,
            const std::vector<double> &posterior_calendar_ages_AD);

private:
    double find_probability_and_ranges_for_cut_off(double cut_off, std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> simplify_ranges(std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> get_ranges_by_intercepts(double probability);
    void _add_text_ranges(const std::vector<double>& values);
};


#endif //CARBONDATE_POSTERIORDENSITYOUTPUT_H
