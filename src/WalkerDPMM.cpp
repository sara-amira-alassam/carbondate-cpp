#include <algorithm>
#include <ctime>
#include <string>
#include <utility>
#include "WalkerDPMM.h"
#include "helpers.h"


void WalkerDPMM::initialise_storage(){
    DPMM::initialise_storage();
    n_clust.resize(n_out);
    weight.resize(n_out);
}

void WalkerDPMM::initialise_clusters() {
    n_weights = n_clust_i = n_clust[0] = std::min(10, n_obs);
    phi_i.resize(n_weights);
    tau_i.resize(n_weights);
    v.resize(n_weights);
    weight_i.resize(n_weights);

    alpha_i = alpha[0] = 2.;
    mu_phi_i = mu_phi[0] = A;

    for (int c = 0; c < n_weights; c++) {
        tau_i[c] = rgamma(nu1, 1./nu2);
    }
    for (int c = 0; c < n_weights; c++) {
        phi_i[c] = rnorm(mu_phi_i, pow(lambda*tau_i[c], -0.5));
    }
    double cumprod = 1.0;
    for (int c = 0; c < n_weights; c++) {
        v[c] = rbeta(1., alpha_i);
        if (c > 0) cumprod *= 1. - v[c - 1];
        weight_i[c] = v[c] * cumprod;
    }

    cluster_ids_i.resize(n_obs);
    get_sample_ids(cluster_ids_i, 1, n_clust_i);
    cluster_ids[0] = cluster_ids_i;
    phi[0] = phi_i;
    tau[0] = tau_i;
    weight[0] = weight_i;
}

void WalkerDPMM::perform_update_step() {
    std::vector<double> u(n_obs); // Auxiliary variables
    double min_u = 1.;            // Minimum value of auxilliary values

    // Create auxiliary variables u
    for (int k = 0; k < n_obs; k++) {
        u[k] = runif(0., weight_i[cluster_ids_i[k] - 1]);
        if (u[k] < min_u) min_u = u[k];
    }

    update_weights(u, min_u);
    update_phi_and_tau();
    update_cluster_ids(u);
    update_n_clust();
    update_alpha();
    update_mu_phi();
    update_calendar_ages();
}

void WalkerDPMM::update_weights(const std::vector<double>& u, double min_u) {
    // Iteratively updates the weights until we have all we need
    // This can change the number of weights
    int cluster_id = 0;
    double brprod = 1., sum_weight = 0.;

    while (sum_weight < 1. - min_u) {
        cluster_id++;
        if (cluster_id <= n_weights) {
            update_v_element(cluster_id, brprod, u);
        } else {
            v.push_back(rbeta(1., alpha_i));   // v_(cluster_id) just from prior beta
            weight_i.push_back(0.); // Reserve space for the new weight
        }
        weight_i[cluster_id - 1] = brprod * v[cluster_id - 1];
        sum_weight += weight_i[cluster_id - 1];
        brprod *= (1. - v[cluster_id - 1]);
    }
    n_weights = cluster_id;
    v.resize(n_weights);
    weight_i.resize(n_weights);
}

void WalkerDPMM::update_v_element(int cluster_id, double brprod, const std::vector<double>& u) {
    std::vector<double> prodv(n_weights);
    bool prodv_set = false;
    double b_sub, b_sub_max = 0.;
    int index;
    double a = 0., b = 1., a_bar, b_bar;

    // We have to work out a and b and sample a new v_(cluster_id) by inverse cdf
    for (int k = 0; k < n_obs; ++k) {
        if (cluster_ids_i[k] > cluster_id) {
            if (!prodv_set) {
                // Vector to aid denominator for b: ith entry is Prod_{l < i, l != j} (1- v_l)
                // In fact strictly here  l == j included, but we multiply by (1 - v_j) below
                prodv[0] = 1.0 - v[0];
                for (int c = 1; c < n_weights; c++) prodv[c] = prodv[c - 1] * (1. - v[c]);
                prodv_set = true;
            }
            index = cluster_ids_i[k] - 1;
            b_sub = u[k] / (v[index] * prodv[index - 1]);
            if (b_sub > b_sub_max) b_sub_max = b_sub;
        } else if ((cluster_ids_i[k] == cluster_id) && (u[k] > a)) {
            // Set 'a' to the max value of u[k] for the observations in the current cluster
            a = u[k];
        }
    }
    a /= brprod;
    b -= b_sub_max * (1.0 - v[cluster_id - 1]);

    a_bar = pow(1. - a, alpha_i);
    b_bar = a_bar - pow(1. - b, alpha_i);

    v[cluster_id - 1] = 1. - pow(a_bar - b_bar * runif(0., 1.), 1. / alpha_i);
}

void WalkerDPMM::update_phi_and_tau() {
    std::vector<double> cluster_calendar_ages(0);
    cluster_calendar_ages.reserve(n_obs);

    phi_i.resize(n_weights);
    tau_i.resize(n_weights);
    for (int c = 1; c <= n_weights; c++) {
        // Find out which observations belong in this cluster
        for (int j = 0; j < n_obs; j++) {
            if (cluster_ids_i[j] == c) cluster_calendar_ages.push_back(calendar_age_i[j]);
        }
        if (cluster_calendar_ages.empty()) {
            // No observations in this cluster, so sample from the prior
            tau_i[c-1] = rgamma(nu1, 1./nu2);
            phi_i[c-1] = rnorm(mu_phi_i, 1./sqrt(lambda * tau_i[c-1]));
        } else {
            // There are some observations, so update the phi and tau for this
            // cluster (conjugate according to NormalGamma prior)
            update_cluster_phi_and_tau(c, cluster_calendar_ages);
            cluster_calendar_ages.clear();
        }
    }
}

void WalkerDPMM::update_cluster_ids(const std::vector<double>& u) {
    std::vector<int> poss_cluster_ids(0);
    poss_cluster_ids.reserve(n_weights);
    std::vector<double> dens(0);
    dens.reserve(n_weights);

    for (int j = 0; j < n_obs; j++) {
        for (int c = 1; c <= n_weights; c++) {
            if (weight_i[c - 1] > u[j]) {
                poss_cluster_ids.push_back(c);
                dens.push_back(dnorm4(phi_i[c - 1], calendar_age_i[j], sqrt(1. / tau_i[c - 1]), 0));
            }
        }
        cluster_ids_i[j] = poss_cluster_ids[sample_integer(poss_cluster_ids.size(), dens, false)];
        poss_cluster_ids.clear();
        dens.clear();
    }
}

void WalkerDPMM::update_n_clust() {
    // Find number of distinct populated clusters
    int cluster_id;
    std::vector<int> observations_per_cluster(n_weights, 0);
    n_clust_i = 0;
    for (int j = 0; j < n_obs; j++) {
        cluster_id = cluster_ids_i[j];
        if (observations_per_cluster[cluster_id - 1] == 0) n_clust_i++;
        observations_per_cluster[cluster_id - 1]++;
    }
}

double WalkerDPMM::alpha_log_likelihood(double alpha_value) {
    return n_clust_i * log(alpha_value) + lgamma(alpha_value) - lgamma(alpha_value + n_obs);
}

double WalkerDPMM::alpha_log_prior(double alpha_value) {
    return dgamma(alpha_value, alpha_shape, 1./alpha_rate, 1);
}

void WalkerDPMM::update_alpha() {
    double updated_alpha = -1.;
    double prop_sd = 1.;        // Standard deviation for sampling proposed value of alpha
    double log_prior_rate, log_likelihood_rate, log_proposal_rate, hr;

    // Sample new alpha from truncated normal distribution
    while (updated_alpha <= 0.) updated_alpha = rnorm(alpha_i, prop_sd);

    log_prior_rate = alpha_log_prior(updated_alpha) - alpha_log_prior(alpha_i);
    log_likelihood_rate = alpha_log_likelihood(updated_alpha) - alpha_log_likelihood(alpha_i);
    // Adjust for non-symmetric truncated normal proposal
    log_proposal_rate = pnorm5(alpha_i, 0., 1., 1, 1) - pnorm5(updated_alpha, 0., 1., 1, 1);
    hr = exp(log_prior_rate + log_likelihood_rate + log_proposal_rate);

    // Accept or reject new alpha
    if (runif(0., 1.) < hr) alpha_i = updated_alpha;
}

void WalkerDPMM::store_current_values(int output_index) {
    DPMM::store_current_values(output_index);
    n_clust[output_index] = n_clust_i;
    weight[output_index] = weight_i;
}

DensityData WalkerDPMM::get_predictive_density(int n_posterior_samples, double resolution, double quantile_edge_width) {
    int n_burn = floor(n_out / 2);
    std::vector<int> sample_ids(n_posterior_samples);
    int s; //current sample id
    get_sample_ids(sample_ids, n_burn - 1, n_out - 1);

    double min_calendar_age = std::numeric_limits<double>::infinity();
    double max_calendar_age = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_out; i++) {
        for (int j = 0; j < n_obs; j++) {
            if (calendar_age[i][j] < min_calendar_age) min_calendar_age = calendar_age[i][j];
            if (calendar_age[i][j] > max_calendar_age) max_calendar_age = calendar_age[i][j];
        }
    }
    min_calendar_age = floor(min_calendar_age);
    max_calendar_age = ceil(max_calendar_age);
    int n_points = (int) ((max_calendar_age - min_calendar_age) / resolution) + 2;
    DensityData density_data(n_points);

    // A vector of vectors to represent the density matrix
    std::vector<std::vector<double>> density_samples(
            n_points, std::vector<double>(n_posterior_samples, 0));
    double sum_weight, logmarg;

    std::vector<double> cal_age_BP(n_points), sum_ages(n_points, 0);
    for (int i = 0; i < n_points; i++) cal_age_BP[i] = min_calendar_age + i * resolution;
    for (int j = 0; j < n_posterior_samples; j++) {
        s = sample_ids[j];
        for (int i = 0; i < n_points; i++) {
            sum_weight = 0.;
            int n_all_clust = (int) weight[s].size();
            for (int c = 0; c < n_all_clust; c++) {
                density_samples[i][j] += weight[s][c]
                                         * dnorm4(cal_age_BP[i], phi[s][c], 1. / sqrt(tau[s][c]), 0);
                sum_weight += weight[s][c];
            }
            // The predictive density for a new observation is a scaled t-distribution
            logmarg = log_marginal_normal_gamma(cal_age_BP[i], mu_phi[s]);
            density_samples[i][j] += (1. - sum_weight) * exp(logmarg);

            sum_ages[i] += density_samples[i][j];
        }
    }

    // Now reverse as we populate the DensityData object as we're converting from calBP to AD
    int i_reversed;
    for (int i = 0; i < n_points; i++) {
        i_reversed = n_points - 1 - i;
        density_data.cal_age_AD[i_reversed] = to_calAD(cal_age_BP[i]);
        density_data.mean[i_reversed] = sum_ages[i] / n_posterior_samples;
        edge_quantiles(
                density_samples[i],
                quantile_edge_width,
                density_data.ci_lower[i_reversed],
                density_data.ci_upper[i_reversed]);
    }

    return density_data;
}
