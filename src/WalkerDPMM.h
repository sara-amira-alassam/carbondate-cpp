/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_WALKERDPMM_H
#define CARBONDATE_WALKERDPMM_H

#include "DPMM.h"

/*
 * See the parent class for the public methods.
 * The private methods defined here implement the Walker method of updating.
 */
class WalkerDPMM : public DPMM {
private:
    // Instant values of DPMM parameters
    // Note that the number of weights can be more than the number of clusters as they may not all be populated
    int n_weights;
    std::vector<double> v, weight_i;

    // Stored values of DPMM parameters
    std::vector<std::vector<double>> weight;

private:
    void _initialise_storage() override;
    void _initialise_clusters() override;
    void _perform_update_step() override;
    void _store_current_values(int i) override;
    void update_weights(const std::vector<double>& u, double min_u);
    void update_v_element(int cluster_id, double brprod, const std::vector<double>& u);
    void update_phi_and_tau();
    void update_cluster_ids(const std::vector<double>& u);
    void update_n_clust();
    double _alpha_log_likelihood(double alpha_value) override;
    double _calculate_density_sample(int sample_id, double calendar_age) override;
};

#endif //CARBONDATE_WALKERDPMM_H
