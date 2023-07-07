#ifndef CARBONDATE_POLYAURNDPMM_H
#define CARBONDATE_POLYAURNDPMM_H
#include "DPMM.h"

#include <utility>

class PolyaUrnDPMM : public DPMM  {
private:
    // Instant values of DPMM parameters
    std::vector<int> observations_per_cluster_i;

    // Stored values of DPMM parameters
    std::vector<std::vector<int>> observations_per_cluster;

private:
    void initialise_storage() override;
    void initialise_clusters() override;
    void perform_update_step() override;
    void store_current_values(int i) override;
    void update_cluster_ids();
    void create_new_phi_and_tau(double calendar_age, double &phi, double &tau);
    void update_phi_and_tau();
    double alpha_log_likelihood(double alpha_value) override;
    double calculate_density_sample(int sample_id, double calendar_age) override;

public:
    explicit PolyaUrnDPMM(std::string file_prefix) : DPMM(std::move(file_prefix)) {};
};


#endif //CARBONDATE_POLYAURNDPMM_H
