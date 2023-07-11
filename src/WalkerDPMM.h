#ifndef CARBONDATE_WALKERDPMM_H
#define CARBONDATE_WALKERDPMM_H
#include "DPMM.h"


class WalkerDPMM : public DPMM {
private:
    // Instant values of DPMM parameters
    // Note that the number of weights can be > number of clusters as may not all be populated
    int n_weights;
    std::vector<double> v, weight_i;

    // Stored values of DPMM parameters
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
    double alpha_log_likelihood(double alpha_value) override;
    double calculate_density_sample(int sample_id, double calendar_age) override;

public:
    explicit WalkerDPMM(std::string file_prefix) : DPMM(std::move(file_prefix)) {};
};

#endif //CARBONDATE_WALKERDPMM_H
