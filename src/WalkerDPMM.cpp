#include <algorithm>
#include <utility>
#include "WalkerDPMM.h"
#include "helpers.h"

void WalkerDPMM::initialise(
        std::vector<double> i_c14_age,
        std::vector<double> i_c14_sig,
        std::vector<double> cc_cal_age,
        std::vector<double> cc_c14_age,
        std::vector<double> cc_c14_sig
) {
    c14_age = std::move(i_c14_age);
    c14_sig = std::move(i_c14_sig);

    n_obs = c14_age.size();
    n_out = 1;

    calcurve.cal_age = std::move(cc_cal_age);
    calcurve.c14_age = std::move(cc_c14_age);
    calcurve.c14_sig = std::move(cc_c14_sig);

    initialise_storage();
    initialise_calendar_age();
    initialise_hyperparameters();
    initialise_clusters();
    interpolate_calibration_curve();
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
    unsigned n_points = calcurve.cal_age.size();
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
    unsigned n = calcurve.cal_age.size();
    std::vector<int> perm(n);
    std::vector<double> sorted_cal_age(calcurve.cal_age.begin(), calcurve.cal_age.end());
    int k_start = 1.; // Start age for interpolated calendar age

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
