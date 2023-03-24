#include <algorithm>
#include <ctime>
#include <string>
#include <utility>
#include "WalkerDPMM.h"
#include "helpers.h"

void WalkerDPMM::initialise(
        std::vector<double> i_c14_age,
        std::vector<double> i_c14_sig,
        std::vector<double> cc_cal_age,
        std::vector<double> cc_c14_age,
        std::vector<double> cc_c14_sig,
        int rng_seed
) {
    if (rng_seed == 0) {
        std::time_t t1, t2;
        t1 = time(nullptr);
        t2 = time(nullptr);
        set_seed(t1, t2);
    } else {
        set_seed(rng_seed, 1);
    }

    c14_age = std::move(i_c14_age);
    c14_sig = std::move(i_c14_sig);

    n_obs = (int) c14_age.size();
    n_out = 1;

    calcurve.cal_age = std::move(cc_cal_age);
    calcurve.c14_age = std::move(cc_c14_age);
    calcurve.c14_sig = std::move(cc_c14_sig);

    interpolate_calibration_curve();
    initialise_storage();
    initialise_calendar_age();
    initialise_hyperparameters();
    initialise_clusters();
}

void WalkerDPMM::initialise_storage(){
    alpha.resize(n_out);
    n_clust.resize(n_out);
    mu_phi.resize(n_out);
    phi.resize(n_out);
    tau.resize(n_out);
    weight.resize(n_out);
    calendar_age.resize(n_out);
    cluster_ids.resize(n_out);
}

void WalkerDPMM::initialise_calendar_age() {
    int n_points = (int) calcurve.cal_age.size();
    double current_prob, max_prob, most_probably_age;
    calendar_age_i.resize(n_obs);

    for (int i = 0; i < n_obs; i++) {
        max_prob = 0.;
        for (int j = 0; j < n_points; j++) {
            current_prob = dnorm4(
                    c14_age[i],
                    calcurve.c14_age[j],
                    sqrt(pow(calcurve.c14_sig[j], 2) + pow(c14_sig[i], 2)),
                    0);
            if (current_prob > max_prob) {
                most_probably_age = calcurve.cal_age[j];
                max_prob = current_prob;
            }
        }
        calendar_age_i[i] = most_probably_age;
    }
    calendar_age[0] = calendar_age_i;
}

void WalkerDPMM::initialise_hyperparameters() {
    double calendar_age_range, calendar_age_prec;

    calendar_age_range = max_diff(calendar_age[0]);
    calendar_age_prec = pow(0.1 * mad(calendar_age[0]), -2);

    A = median(calendar_age[0]);
    B = 1./pow(calendar_age_range, 2);
    mu_phi[0] = A;

    lambda = pow(100. / calendar_age_range, 2);
    nu1 = 0.25;
    nu2 = nu1 / calendar_age_prec;

    slice_width = std::max(1000., max_diff(c14_age) / 2.);
}

void WalkerDPMM::initialise_clusters() {
    n_weights = n_clust_i = n_clust[0] = 10;
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
    get_sample_ids(cluster_ids_i, 1, n_clust_i, n_obs);
    cluster_ids[0] = cluster_ids_i;
    phi[0] = phi_i;
    tau[0] = tau_i;
    weight[0] = weight_i;
}

void WalkerDPMM::interpolate_calibration_curve() {
   int n = (int) calcurve.cal_age.size();
    std::vector<int> perm(n);
    std::vector<double> sorted_cal_age(calcurve.cal_age.begin(), calcurve.cal_age.end());
    int k_start = 1.; // TODO : allow different start ages
    // Start age for interpolated calendar age

    yearwise_calcurve.cal_age.resize(max_year_bp);
    yearwise_calcurve.c14_age.resize(max_year_bp);
    yearwise_calcurve.c14_sig.resize(max_year_bp);

    for (int i = 0; i < n; i++) perm[i] = i;
    rsort_with_index(&sorted_cal_age[0], &perm[0], n);

    int i = 0;
    while (sorted_cal_age[i] < k_start) {
        i++;
    }
    for (int k = k_start; k <= max_year_bp; k++) {
        yearwise_calcurve.cal_age[k - 1] = k;
        if (sorted_cal_age[i] == k) {
            yearwise_calcurve.c14_age[k - 1] = calcurve.c14_age[perm[i]];
            yearwise_calcurve.c14_sig[k - 1] = calcurve.c14_sig[perm[i]];
        } else if (sorted_cal_age[i + 1] == k) {
            i++;
            yearwise_calcurve.c14_age[k - 1] = calcurve.c14_age[perm[i]];
            yearwise_calcurve.c14_sig[k - 1] = calcurve.c14_sig[perm[i]];
        } else {
            yearwise_calcurve.c14_age[k - 1] = interpolate_linear(
                    k,
                    sorted_cal_age[i],
                    sorted_cal_age[i + 1],
                    calcurve.c14_age[perm[i]],
                    calcurve.c14_age[perm[i + 1]]
            );
            yearwise_calcurve.c14_sig[k - 1] = interpolate_linear(
                    k,
                    sorted_cal_age[i],
                    sorted_cal_age[i + 1],
                    calcurve.c14_sig[perm[i]],
                    calcurve.c14_sig[perm[i + 1]]
            );
        }
    }
}

void WalkerDPMM::calibrate(int n_iter, int n_thin) {
    n_out = n_iter/n_thin + 1;
    initialise_storage();
    for (int i = 1; i <= n_iter; i++) {
        perform_update_step();
        if (i % n_thin == 0) {
            update_progress_bar(i * 1. / n_iter);
            store_current_values(i / n_thin);
        }
    }
    printf("\n");
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
    std::vector<double> prodv(n_weights + 1);   // TODO: does this need to be +1??
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
            cluster_calendar_ages.resize(0);
        }
    }
}

void WalkerDPMM::update_cluster_phi_and_tau(
        int cluster_id, const std::vector<double>& cluster_calendar_ages) {
    int num_in_cluster = (int) cluster_calendar_ages.size();
    std::vector<double> calendar_age_diff(num_in_cluster);
    double calendar_age_mean;
    double s;
    double mu_phi_new;
    double lambda_new, nu1_new, nu2_new;

    calendar_age_mean = mean(cluster_calendar_ages);
    for (unsigned i = 0; i < num_in_cluster; ++i) {
        calendar_age_diff[i] = pow(cluster_calendar_ages[i] - calendar_age_mean, 2);
    }
    s = mean(calendar_age_diff);

    // Update parameters according to conjugate prior
    nu1_new = nu1 + num_in_cluster / 2.;
    nu2_new = nu2;
    nu2_new += 0.5 * num_in_cluster
               * (s + lambda * pow(calendar_age_mean - mu_phi_i, 2)/(lambda + num_in_cluster));

    lambda_new = lambda + num_in_cluster;
    mu_phi_new = (lambda * mu_phi_i + num_in_cluster * calendar_age_mean) / (lambda + num_in_cluster);

    tau_i[cluster_id - 1] = rgamma(nu1_new, 1./nu2_new);
    phi_i[cluster_id - 1] = rnorm(mu_phi_new, 1./sqrt(lambda_new * tau_i[cluster_id - 1]));
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
        poss_cluster_ids.resize(0);
        dens.resize(0);
    }
}

void WalkerDPMM::update_n_clust() {
    // Find number of distinct populated clusters
    int cluster_id;
    std::vector<int> observations_per_cluster(2*n_obs, 0);  // Allocate plenty of space
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

void WalkerDPMM::update_mu_phi() {
    // Updates mu_phi via Gibbs sampling based upon current (phi, tau) values
    double posterior_mean, posterior_precision;
    double sum_tau = 0.;
    double sum_tau_mult_phi = 0.;

    for (int c = 1; c <= n_weights; c++) {
        sum_tau += tau_i[c - 1];
        sum_tau_mult_phi += tau_i[c - 1] * phi_i[c - 1];
    }
    posterior_mean = (A * B + lambda * sum_tau_mult_phi) / (B + lambda * sum_tau);
    posterior_precision = B + lambda * sum_tau;
    mu_phi_i = rnorm(posterior_mean, 1. / sqrt(posterior_precision));
}

double WalkerDPMM::cal_age_log_likelihood(
        double cal_age,
        double cluster_mean,
        double cluster_sig,
        double obs_c14_age,
        double obs_c14_sig
) {
    double loglik;
    double cc_c14_age, cc_c14_sig;
    int yr;

    yr = (int) cal_age - 1;
    // TODO: Change this
    if ((yr < 0) || (yr >= max_year_bp)) {  // out of range
        return -std::numeric_limits<double>::infinity();
    }
    cc_c14_age = yearwise_calcurve.c14_age[yr];
    cc_c14_sig= yearwise_calcurve.c14_sig[yr];

    loglik = dnorm4(cal_age, cluster_mean, cluster_sig, 1);
    loglik += dnorm4(
            obs_c14_age, cc_c14_age, sqrt(cc_c14_sig*cc_c14_sig + obs_c14_sig*obs_c14_sig), 1);

    return loglik;
}

void WalkerDPMM::update_calendar_ages() {
    // updates calendar ages using slice sampling
    int cluster_id;
    double cluster_mean, cluster_sig;
    double y;     // Slice height
    double L, R;  // Each side of calendar age slice interval
    double J, K;  // Max steps on left and right hand sides
    double x0, x1;    // Old value and current sampled value of calendar age

    for (int j = 0; j < n_obs; j++) {
        cluster_id = cluster_ids_i[j];
        cluster_mean = phi_i[cluster_id - 1];
        cluster_sig = 1.0 / sqrt(tau_i[cluster_id - 1]);
        x0 = calendar_age_i[j];

        // Slice height
        y = cal_age_log_likelihood(x0, cluster_mean, cluster_sig, c14_age[j], c14_sig[j]) - rexp(1);

        //////////////////////////////////////////////
        // Find the slice interval

        // Initialise slice
        L = x0 - slice_width * runif(0., 1.);
        R = L + slice_width;

        // Decide how far we can extend on both sides
        J = floor(slice_multiplier * runif(0., 1.));
        K = slice_multiplier - 1. - J;

        // LHS stepping out
        while ((J > 0) &&
               (y < cal_age_log_likelihood(L, cluster_mean, cluster_sig, c14_age[j], c14_sig[j]))) {
            L -= slice_width;
            J -= 1.;
        }

        // RHS stepping out
        while ((K > 0) &&
               (y < cal_age_log_likelihood(R, cluster_mean, cluster_sig, c14_age[j], c14_sig[j]))) {
            R += slice_width;
            K -= 1.;
        }

        //////////////////////////////////////////////
        // Get sampled value from the slice interval
        while (true) {
            x1 = L + runif(0., 1.) * (R - L);

            // Break loop if we have sampled satisfactorily
            if (y < cal_age_log_likelihood(x1, cluster_mean, cluster_sig, c14_age[j], c14_sig[j])) {
                calendar_age_i[j] = x1;
                break;
            }

            // Else shrink the interval
            if (x1 < x0) {
                L = x1;
            } else {
                R = x1;
            }
        }
    }
}

void WalkerDPMM::store_current_values(int output_index) {
    alpha[output_index] = alpha_i;
    mu_phi[output_index] = mu_phi_i;
    n_clust[output_index] = n_clust_i;
    phi[output_index] = phi_i;
    tau[output_index] = tau_i;
    weight[output_index] = weight_i;
    calendar_age[output_index] = calendar_age_i;
    cluster_ids[output_index] = cluster_ids_i;
}

// Function which works out the marginal of a calendar age when
// theta ~ N(phi, sd = sqrt(1/tau)) and (phi,tau) are NormalGamma
double WalkerDPMM::log_marginal_normal_gamma(double cal_age, double mu_phi_s) {
    double logden, margprec, margdf;

    margprec = (nu1 * lambda) / (nu2 * (lambda + 1.));
    margdf = 2. * nu1;

    logden = lgamma((margdf + 1.) / 2.) - lgamma(margdf / 2.);
    logden += 0.5 * (log(margprec) - log(margdf) - log(M_PI));
    logden -= ((margdf + 1) / 2) * log(1 + margprec * pow(cal_age - mu_phi_s, 2) / margdf);

    return logden;
}

DensityData WalkerDPMM::get_predictive_density(
        int n_posterior_samples,
        double resolution,
        double quantile_edge_width) {

    int n_burn = floor(n_out / 2);
    std::vector<int> sample_ids(n_posterior_samples);
    int s; //current sample id
    get_sample_ids(sample_ids, n_burn - 1, n_out - 1, n_posterior_samples);

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

std::vector<double> WalkerDPMM::get_posterior_calendar_ages(int ident) {
    int n_burn = n_out / 2;
    int n_count = n_out - n_burn;

    std::vector<double> posterior_calendar_ages(n_count);
    for (int i = 0; i < n_count; i++) {
        posterior_calendar_ages[i] = to_calAD(calendar_age[i + n_burn][ident]);
    }

    return posterior_calendar_ages;
}
