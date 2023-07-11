#include "DPMM.h"
#include "helpers.h"
#include "work.h"

void DPMM::initialise(
        std::vector<double> i_rc_determinations,
        std::vector<double> i_rc_sigmas,
        bool i_f14c_inputs,
        std::vector<double> cc_cal_age,
        std::vector<double> cc_c14_age,
        std::vector<double> cc_c14_sig,
        int rng_seed) {
    if (rng_seed == 0) {
        std::time_t t1, t2;
        t1 = time(nullptr);
        t2 = time(nullptr);
        set_seed(t1, t2);
    } else {
        set_seed(rng_seed, 1);
    }

    rc_determinations = std::move(i_rc_determinations);
    rc_sigmas = std::move(i_rc_sigmas);
    f14c_inputs = i_f14c_inputs;

    n_obs = (int) rc_determinations.size();
    n_out = 1;

    calcurve.cal_age = std::move(cc_cal_age);
    calcurve.c14_age = std::move(cc_c14_age);
    calcurve.c14_sig = std::move(cc_c14_sig);
    convert_to_f14c_age(calcurve.c14_age, calcurve.c14_sig, calcurve.f14c_age, calcurve.f14c_sig);

    interpolate_calibration_curve();
    initialise_storage();
    initialise_calendar_age_and_spd_ranges();
    initialise_hyperparameters();
    initialise_clusters();
}

void DPMM::initialise_storage() {
    calendar_age.resize(n_out);
    alpha.resize(n_out);
    mu_phi.resize(n_out);
    phi.resize(n_out);
    tau.resize(n_out);
    n_clust.resize(n_out);
    cluster_ids.resize(n_out);
}

void DPMM::initialise_calendar_age_and_spd_ranges() {
    int n_points = (int) yearwise_calcurve.cal_age.size();
    double sum_prob, max_prob, most_probably_age;
    std::vector<double> current_prob(n_points), spd(n_points), cum_prob_spd(n_points);
    calendar_age_i.resize(n_obs);

    for (int i = 0; i < n_obs; i++) {
        sum_prob = max_prob = 0.;
        for (int j = 0; j < n_points; j++) {
            current_prob[j] = dnorm4(
                    rc_determinations[i],
                    yearwise_calcurve.rc_age[j],
                    sqrt(pow(yearwise_calcurve.rc_sig[j], 2) + pow(rc_sigmas[i], 2)),
                    0);
            sum_prob += current_prob[j];
            if (current_prob[j] > max_prob) {
                most_probably_age = yearwise_calcurve.cal_age[j];
                max_prob = current_prob[j];
            }
        }
        calendar_age_i[i] = most_probably_age;
        for (int j = 0; j < n_points; j++) {
            spd[j] += current_prob[j] / sum_prob;
        }
    }
    calendar_age[0] = calendar_age_i;

    cum_prob_spd[0] = spd[0] / n_obs;
    for (int i = 1; i < n_points; i++) cum_prob_spd[i] = cum_prob_spd[i - 1] + spd[i] / n_obs;

    spd_range_1_sigma[0] = yearwise_calcurve.cal_age[get_left_boundary(cum_prob_spd, (1 - 0.683)/2.)];
    spd_range_1_sigma[1] = yearwise_calcurve.cal_age[get_right_boundary(cum_prob_spd, (1 + 0.683)/2.)];
    spd_range_2_sigma[0] = yearwise_calcurve.cal_age[get_left_boundary(cum_prob_spd, (1 - 0.954)/2.)];
    spd_range_2_sigma[1] = yearwise_calcurve.cal_age[get_right_boundary(cum_prob_spd, (1 + 0.954)/2.)];
    spd_range_3_sigma[0] = yearwise_calcurve.cal_age[get_left_boundary(cum_prob_spd, (1 - 0.997)/2.)];
    spd_range_3_sigma[1] = yearwise_calcurve.cal_age[get_right_boundary(cum_prob_spd, (1 + 0.997)/2.)];
}

void DPMM::initialise_hyperparameters() {
    double calendar_age_range, calendar_age_prec;

    A = mean(spd_range_2_sigma);
    B = 1./pow(max_diff(spd_range_2_sigma), 2);

    calendar_age_range = 0.05 * max_diff(spd_range_1_sigma);
    calendar_age_prec = pow(calendar_age_range, -2);

    lambda = pow(100. / max_diff(spd_range_3_sigma), 2);
    nu1 = 0.25;
    nu2 = nu1 / calendar_age_prec;

    slice_width = max_diff(spd_range_3_sigma);
}

void DPMM::initialise_clusters() {}

void DPMM::interpolate_calibration_curve() {
    int n = (int) calcurve.cal_age.size();
    std::vector<int> perm(n);
    std::vector<double> sorted_cal_age(calcurve.cal_age.begin(), calcurve.cal_age.end());
    std::vector<double> *calcurve_rc_age, *calcurve_rc_sig;

    for (int i = 0; i < n; i++) perm[i] = i;
    rsort_with_index(&sorted_cal_age[0], &perm[0], n);

    // Start and end date for interpolated calendar age
    int y_start = (int) floor(sorted_cal_age[0]);
    int y_end = (int) ceil(sorted_cal_age[n - 1]);
    int n_interp = y_end - y_start + 1;

    yearwise_calcurve.cal_age.resize(n_interp);
    yearwise_calcurve.rc_age.resize(n_interp);
    yearwise_calcurve.rc_sig.resize(n_interp);

    if (f14c_inputs) {
        calcurve_rc_age = &calcurve.f14c_age;
        calcurve_rc_sig = &calcurve.f14c_sig;
    } else {
        calcurve_rc_age = &calcurve.c14_age;
        calcurve_rc_sig = &calcurve.c14_sig;
    }

    int i = 0;
    for (int k = 0; k < n_interp; k++) {
        yearwise_calcurve.cal_age[k] = k + y_start;
        if (sorted_cal_age[i] == yearwise_calcurve.cal_age[k]) {
            yearwise_calcurve.rc_age[k] = calcurve_rc_age->at(perm[i]);
            yearwise_calcurve.rc_sig[k] = calcurve_rc_sig->at(perm[i]);
        } else if (sorted_cal_age[i + 1] == yearwise_calcurve.cal_age[k]) {
            i++;
            yearwise_calcurve.rc_age[k] = calcurve_rc_age->at(perm[i]);
            yearwise_calcurve.rc_sig[k] = calcurve_rc_sig->at(perm[i]);
        } else {
            yearwise_calcurve.rc_age[k] = interpolate_linear(
                    yearwise_calcurve.cal_age[k],
                    sorted_cal_age[i],
                    sorted_cal_age[i + 1],
                    calcurve_rc_age->at(perm[i]),
                    calcurve_rc_age->at(perm[i + 1])
            );
            yearwise_calcurve.rc_sig[k] = interpolate_linear(
                    yearwise_calcurve.cal_age[k],
                    sorted_cal_age[i],
                    sorted_cal_age[i + 1],
                    calcurve_rc_sig->at(perm[i]),
                    calcurve_rc_sig->at(perm[i + 1])
            );
        }
    }
}

void DPMM::calibrate(int n_iter, int n_thin) {
    n_out = n_iter/n_thin + 1;
    initialise_storage();
    for (int i = 1; i <= n_iter; i++) {
        check_for_work_file(_file_prefix);
        perform_update_step();
        if (i % n_thin == 0) {
            update_progress_bar(i * 1. / n_iter);
            store_current_values(i / n_thin);
        }
        if (i % n_work_update == 0) {
            update_work_file_mcmc(_file_prefix, double (i) / n_iter, i);
        }
    }
}

void DPMM::store_current_values(int output_index) {
    calendar_age[output_index] = calendar_age_i;
    alpha[output_index] = alpha_i;
    mu_phi[output_index] = mu_phi_i;
    phi[output_index] = phi_i;
    tau[output_index] = tau_i;
    cluster_ids[output_index] = cluster_ids_i;
    n_clust[output_index] = n_clust_i;
}

double DPMM::cal_age_log_likelihood(
        double cal_age,
        double cluster_mean,
        double cluster_sig,
        double obs_c14_age,
        double obs_c14_sig
) {
    double loglik;
    double cc_c14_age, cc_c14_sig;
    int yr_index = (int) (cal_age - yearwise_calcurve.cal_age[0]);

    if ((yr_index < 0) || (yr_index >= yearwise_calcurve.rc_age.size())) {  // out of range
        return -std::numeric_limits<double>::infinity();
    }
    cc_c14_age = yearwise_calcurve.rc_age[yr_index];
    cc_c14_sig= yearwise_calcurve.rc_sig[yr_index];

    loglik = dnorm4(cal_age, cluster_mean, cluster_sig, 1);
    loglik += dnorm4(
            obs_c14_age, cc_c14_age, sqrt(cc_c14_sig*cc_c14_sig + obs_c14_sig*obs_c14_sig), 1);

    return loglik;
}

void DPMM::update_calendar_ages() {
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
        y = cal_age_log_likelihood(x0, cluster_mean, cluster_sig, rc_determinations[j], rc_sigmas[j]) - rexp(1);

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
               (y < cal_age_log_likelihood(L, cluster_mean, cluster_sig, rc_determinations[j], rc_sigmas[j]))) {
            L -= slice_width;
            J -= 1.;
        }

        // RHS stepping out
        while ((K > 0) &&
               (y < cal_age_log_likelihood(R, cluster_mean, cluster_sig, rc_determinations[j], rc_sigmas[j]))) {
            R += slice_width;
            K -= 1.;
        }

        //////////////////////////////////////////////
        // Get sampled value from the slice interval
        while (true) {
            x1 = L + runif(0., 1.) * (R - L);

            // Break loop if we have sampled satisfactorily
            if (y < cal_age_log_likelihood(x1, cluster_mean, cluster_sig, rc_determinations[j], rc_sigmas[j])) {
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

// Function which works out the marginal of a calendar age when
// theta ~ N(phi, sd = sqrt(1/tau)) and (phi, tau) are NormalGamma
double DPMM::log_marginal_normal_gamma(double cal_age, double mu_phi_s) {
    double logden, margprec, margdf;

    margprec = (nu1 * lambda) / (nu2 * (lambda + 1.));
    margdf = 2. * nu1;

    logden = lgamma((margdf + 1.) / 2.) - lgamma(margdf / 2.);
    logden += 0.5 * (log(margprec) - log(margdf) - log(M_PI));
    logden -= ((margdf + 1) / 2) * log(1 + margprec * pow(cal_age - mu_phi_s, 2) / margdf);

    return logden;
}

void DPMM::update_cluster_phi_and_tau(int cluster_id, const std::vector<double>& cluster_calendar_ages) {
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

void DPMM::update_mu_phi() {
    // Updates mu_phi via Gibbs sampling based upon current (phi, tau) values
    double posterior_mean, posterior_precision;
    double sum_tau = 0.;
    double sum_tau_mult_phi = 0.;

    for (int c = 1; c <= phi_i.size(); c++) {
        sum_tau += tau_i[c - 1];
        sum_tau_mult_phi += tau_i[c - 1] * phi_i[c - 1];
    }
    posterior_mean = (A * B + lambda * sum_tau_mult_phi) / (B + lambda * sum_tau);
    posterior_precision = B + lambda * sum_tau;
    mu_phi_i = rnorm(posterior_mean, 1. / sqrt(posterior_precision));
}

void DPMM::update_alpha() {
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

double DPMM::alpha_log_prior(double alpha_value) {
    return dgamma(alpha_value, alpha_shape, 1./alpha_rate, 1);
}

double DPMM::alpha_log_likelihood(double alpha_value) {
    // Must be implemented in child class
    return 0.;
}


DensityData DPMM::get_predictive_density(int n_posterior_samples, double resolution, double quantile_edge_width) {
    int n_burn = floor(n_out / 2);
    std::vector<int> sample_ids(n_posterior_samples);
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
    std::vector<std::vector<double>> density_samples(n_points, std::vector<double>(n_posterior_samples, 0));

    std::vector<double> cal_age_BP(n_points), sum_ages(n_points, 0);
    for (int i = 0; i < n_points; i++) cal_age_BP[i] = min_calendar_age + i * resolution;
    for (int j = 0; j < n_posterior_samples; j++) {
        for (int i = 0; i < n_points; i++) {
            density_samples[i][j] = calculate_density_sample(sample_ids[j], cal_age_BP[i]);
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

double DPMM::calculate_density_sample(int sample_id, double calendar_age_BP) {
    // Must be implemented in child class
    return 0.;
}

std::vector<double> DPMM::get_posterior_calendar_ages(int ident) {
    int n_burn = n_out / 2;
    int n_count = n_out - n_burn;

    std::vector<double> posterior_calendar_ages(n_count);
    for (int i = 0; i < n_count; i++) {
        posterior_calendar_ages[i] = to_calAD(calendar_age[i + n_burn][ident]);
    }

    return posterior_calendar_ages;
}
