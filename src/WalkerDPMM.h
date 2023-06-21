#ifndef CARBONDATE_WALKERDPMM_H
#define CARBONDATE_WALKERDPMM_H
#include "DPMM.h"


class WalkerDPMM : public DPMM {
private:
    // Instant values of DPMM parameters and output
    // Note that the number of weights can be > number of clusters as may not all be populated
    // Here n_clust refers to the number of unique cluster ids for the observations
    int n_clust_i, n_weights;
    std::vector<double> v, weight_i;

    // Stored values of DPMM parameters and output
    std::vector<int> n_clust;
    std::vector<std::vector<double>> weight;

private:
    void initialise_storage() override;
    void initialise_clusters() override;
    void perform_update_step() override;
    void store_current_values(int i) override;
    void update_weights(const std::vector<double>& u, double min_u);
    void update_v_element(int cluster_id, double brprod, const std::vector<double>& u);
    void update_phi_and_tau();
    void update_cluster_ids(const std::vector<double>& u);
    void update_n_clust();
    void update_alpha();
    double alpha_log_prior(double alpha_value);
    double alpha_log_likelihood(double alpha_value);

public:
    DensityData get_predictive_density(
            int n_posterior_samples, double resolution, double quantile_edge_width) override;

};

#endif //CARBONDATE_WALKERDPMM_H
