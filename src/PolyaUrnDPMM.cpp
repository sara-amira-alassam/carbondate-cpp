#include "PolyaUrnDPMM.h"

void PolyaUrnDPMM::initialise_storage(){
    DPMM::initialise_storage();
    observations_per_cluster.resize(n_out);
}

void PolyaUrnDPMM::initialise_clusters() {
    n_clust_i = n_clust[0] = std::min(10, n_obs);
    phi_i.resize(n_clust_i);
    tau_i.resize(n_clust_i);

    alpha_i = alpha[0] = 0.0001;
    mu_phi_i = mu_phi[0] = A;

    for (int c = 0; c < n_clust_i; c++) {
        phi_i[c] = rnorm(mu_phi_i, 0.5 * max_diff(spd_range_3_sigma));
        tau_i[c] = pow(max_diff(spd_range_3_sigma) / 4., -2);
    }
    phi[0] = phi_i;
    tau[0] = tau_i;

    cluster_ids_i.resize(n_obs);
    observations_per_cluster_i.resize(n_clust_i);
    int n_populated_clusters = 0;
    while (n_populated_clusters != n_clust_i) {
        for (int & cluster : observations_per_cluster_i) cluster = 0;
        get_sample_ids(cluster_ids_i, 1, n_clust_i);

        for (int id : cluster_ids_i ) {
            if (observations_per_cluster_i[id - 1] == 0) n_populated_clusters++;
            observations_per_cluster_i[id - 1]++;
        }
    }
    cluster_ids[0] = cluster_ids_i;
    observations_per_cluster[0] = observations_per_cluster_i;
}

void PolyaUrnDPMM::perform_update_step() {
    update_cluster_ids();
    update_phi_and_tau();
    update_mu_phi();
    update_calendar_ages();
    update_alpha();
}

void PolyaUrnDPMM::store_current_values(int output_index) {
    DPMM::store_current_values(output_index);
    observations_per_cluster[output_index] = observations_per_cluster_i;
}

void PolyaUrnDPMM::update_cluster_ids() {
    int cluster_id, new_cluster_id;
    std::vector<double> cluster_prob(n_clust_i + 1); // Probability of sampling each cluster ID
    double logmarg;
    double phi_new, tau_new;           // New phi and tau values when we sample a new cluster
    int n_non_empty_clust = n_clust_i;   // Keeps track of how many clusters we have each iteration
    std::vector<int> cluster_id_map;   // Maps the old cluster ids to the new shifted cluster ids

    // Reserve memory to increase vectors if we introduce new clusters
    cluster_prob.reserve(2 * n_clust_i);
    phi_i.reserve(2 * n_clust_i);
    tau_i.reserve(2 * n_clust_i);

    for (int i = 0; i < n_obs; i++) {
        cluster_id = cluster_ids_i[i];

        // Updating this cluster ID, so remove from count
        observations_per_cluster_i[cluster_id - 1]--;
        if (observations_per_cluster_i[cluster_id - 1] == 0) n_non_empty_clust--;

        // Calculate the probability for sampling each new cluster_id for this observation
        // (including a new cluster ID)
        for (int c = 1; c <= n_clust_i; c++) {
            if (observations_per_cluster_i[c - 1] == 0) {
                cluster_prob[c - 1] = 0.;
            } else {
                cluster_prob[c - 1] = dnorm4(calendar_age_i[i], phi_i[c - 1], 1./sqrt(tau_i[c - 1]), 0);
                cluster_prob[c - 1] *= observations_per_cluster_i[c - 1];
            }
            logmarg = log_marginal_normal_gamma(calendar_age_i[i], mu_phi_i);
            cluster_prob[n_clust_i] = exp(logmarg) * alpha_i;
        }
        // Sample cluster ID for the new cluster
        new_cluster_id = sample_integer(n_clust_i + 1, cluster_prob, true);

        // If we've sampled a new cluster, add new phi and tau and update other variables to reflect
        if (new_cluster_id == n_clust_i + 1) {
            create_new_phi_and_tau(calendar_age_i[i], phi_new, tau_new);
            phi_i.push_back(phi_new);
            tau_i.push_back(tau_new);
            cluster_prob.push_back(0.);
            observations_per_cluster_i.push_back(1);
            n_clust_i++;
            n_non_empty_clust++;
        } else {
            observations_per_cluster_i[new_cluster_id - 1]++;
        }
        cluster_ids_i[i] = new_cluster_id;
    }

    // Create a map of old cluster labelling to new cluster labelling
    // Shift phi and tau to remove values where there are no observations
    cluster_id_map.resize(n_clust_i + 1);
    int new_c = 1;
    for (int c = 1; c <= n_clust_i; c++) {
        if (observations_per_cluster_i[c - 1] > 0) {
            cluster_id_map[c] = new_c;
            phi_i[new_c - 1] = phi_i[c - 1];
            tau_i[new_c - 1] = tau_i[c - 1];
            observations_per_cluster_i[new_c - 1] = observations_per_cluster_i[c - 1];
            new_c++;
        }
    }
    phi_i.resize(n_non_empty_clust);
    tau_i.resize(n_non_empty_clust);
    observations_per_cluster_i.resize(n_non_empty_clust);
    n_clust_i = n_non_empty_clust;

    // Change the cluster ID labelling so that there are no skipped values
    for (int i = 0; i < n_obs; i++) cluster_ids_i[i] = cluster_id_map[cluster_ids_i[i]];
}


void PolyaUrnDPMM::create_new_phi_and_tau(double calendar_age, double &phi, double &tau) {

    double nu1_l = nu1 + 0.5;
    double nu2_l = nu2 + lambda * pow(calendar_age - mu_phi_i, 2) / (2. * (lambda + 1.));
    double mu_phi_l = (lambda * mu_phi_i + calendar_age) / (lambda + 1);
    double lambda_l = lambda + 1;

    tau = rgamma(nu1_l, 1. / nu2_l);
    phi = rnorm(mu_phi_l, 1. / sqrt(lambda_l * tau));
}

void PolyaUrnDPMM::update_phi_and_tau() {
    std::vector<double> cluster_calendar_ages(0);
    cluster_calendar_ages.reserve(n_obs);

    for (int c = 1; c <= n_clust_i; c++) {
        // Find out which observations belong in this cluster
        for (int j = 0; j < n_obs; j++) if (cluster_ids_i[j] == c) cluster_calendar_ages.push_back(calendar_age_i[j]);
        update_cluster_phi_and_tau(c, cluster_calendar_ages);
        cluster_calendar_ages.clear();
    }
}

double PolyaUrnDPMM::alpha_log_likelihood(double alpha_value){
    double log_likelihood = n_clust_i * log(alpha_value);
    int c_shift;

    for (int i = 0; i < n_clust_i; i++) {
        // Use max of value or 1 since 0! is 1
        c_shift = observations_per_cluster_i[i] <= 1 ? 1 : observations_per_cluster_i[i] - 1;
        for (int j = 1; j <= c_shift; j++) log_likelihood += log(j);
    }
    for (int i = 0; i < n_obs; i++) log_likelihood -= log(alpha_value + i);
    return log_likelihood;
}

double PolyaUrnDPMM::calculate_density_sample(int sample_id, double calendar_age_BP) {
    double logmarg, density_sample = 0.;

    for (int c = 0; c < n_clust[sample_id]; c++) {
        density_sample += observations_per_cluster[sample_id][c]
                                 * dnorm4(calendar_age_BP, phi[sample_id][c], 1. / sqrt(tau[sample_id][c]), 0);
    }
    // The predictive density for a new observation is a scaled t-distribution
    logmarg = log_marginal_normal_gamma(calendar_age_BP, mu_phi[sample_id]);
    density_sample += alpha[sample_id] * exp(logmarg);
    density_sample /= n_obs + alpha[sample_id];

    return density_sample;
}
