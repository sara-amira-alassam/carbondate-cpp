/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_PREDICTIVEDENSITYOUTPUT_H
#define CARBONDATE_PREDICTIVEDENSITYOUTPUT_H

#include "DensityOutput.h"

class PredictiveDensityOutput : public DensityOutput {
    std::vector<double> _ci_lower;
    std::vector<double> _ci_upper;
    int _n_obs;
    std::string _name;

    void _generate_output_lines() override;
    void set_confidence_intervals(const std::vector<double>& ci_lower, const std::vector<double>& ci_upper);
    std::string _range_lines(int &comment_index) override;

public:
    PredictiveDensityOutput(
            int n_obs,
            int offset,
            double resolution,
            std::string model_name,
            const std::vector<double>& cal_age_AD,
            const std::vector<double>& mean,
            const std::vector<double>& ci_lower,
            const std::vector<double>& ci_upper);
};


#endif //CARBONDATE_PREDICTIVEDENSITYOUTPUT_H
