/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_POLYAURNDPMM_H
#define CARBONDATE_POLYAURNDPMM_H
#include "DPMM.h"

/*
 * See the parent class for the public methods.
 * The private methods defined here implement the Polya Urn method of updating.
 */
class PolyaUrnDPMM : public DPMM  {
private:
    // Instant values of DPMM parameters
    std::vector<int> observations_per_cluster_i;

    // Stored values of DPMM parameters
    std::vector<std::vector<int>> observations_per_cluster;

private:
    void _initialise_storage() override;
    void _initialise_clusters() override;
    void _perform_update_step() override;
    void _store_current_values(int i) override;
    void update_cluster_ids();
    void create_new_phi_and_tau(double calendar_age, double &phi, double &tau);
    void update_phi_and_tau();
    double _alpha_log_likelihood(double alpha_value) override;
    double _calculate_density_sample(int sample_id, double calendar_age) override;
};

#endif //CARBONDATE_POLYAURNDPMM_H
