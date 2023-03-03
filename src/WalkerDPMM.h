#ifndef CARBONDATE_WALKERDPMM_H
#define CARBONDATE_WALKERDPMM_H
#include <vector>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "DensityOutput.h"

struct CalCurve {
    std::vector<double> cal_age;
    std::vector<double> c14_age;
    std::vector<double> c14_sig;
};

struct DensityData {
    std::vector<double> cal_age;
    std::vector<double> mean;
    std::vector<double> ci_lower;
    std::vector<double> ci_upper;
};

class WalkerDPMM {
    std::vector<double> c14_age;  // observed c14 determinations
    std::vector<double> c14_sig;  // c14 determination uncertainties
    CalCurve calcurve;            // original calibration curve data
    CalCurve yearwise_calcurve;   // calibration curve interpolated for every year of calendar age
    int max_year_bp = 50000;      // maximum year for interpolated calendar age

    int n_obs, n_out;

    // Hyperparameters
    double lambda, nu1, nu2;
    double A, B;
    double alpha_shape = 1., alpha_rate = 1.;
    double slice_width, slice_multiplier = 10.;

    // Instant values of DPMM parameters and output
    // Note that the number of weights can be > number of clusters as may not all be populated
    // Here n_clust refers to the number of unique cluster ids for the observations
    double alpha_i, mu_phi_i;
    int n_clust_i, n_weights;
    std::vector<double> phi_i, tau_i, v, weight_i, calendar_age_i;
    std::vector<int> cluster_ids_i;

    // Stored values of DPMM parameters and output
    std::vector<double> alpha, mu_phi;
    std::vector<int> n_clust;
    std::vector<std::vector<double>> phi, tau, weight, calendar_age;
    std::vector<std::vector<int>> cluster_ids;

private:
    void initialise_storage();
    void initialise_calendar_age();
    void initialise_hyperparameters();
    void initialise_clusters();
    void interpolate_calibration_curve();
    void perform_update_step();
    void store_current_values(int i);
    void update_weights(const std::vector<double>& u, double min_u);
    void update_v_element(int cluster_id, double brprod, const std::vector<double>& u);
    void update_phi_and_tau();
    void update_cluster_phi_and_tau(int cluster_id, const std::vector<double>& cluster_calendar_ages);
    void update_cluster_ids(const std::vector<double>& u);
    void update_n_clust();
    void update_alpha();
    void update_mu_phi();
    void update_calendar_ages();
    double alpha_log_prior(double alpha_value);
    double alpha_log_likelihood(double alpha_value);
    double cal_age_log_likelihood(
            double cal_age,
            double prmean,
            double prsig,
            double obs_c14_age,
            double obs_c14_sig
    );
    double log_marginal_normal_gamma(double cal_age, double mu_phi_s);

public:
    void initialise(
            std::vector<double> i_c14_age,
            std::vector<double> i_c14_sig,
            std::vector<double> cc_cal_age,
            std::vector<double> cc_c14_age,
            std::vector<double> cc_c14_sig,
            int rng_seed = 0
    );
    void calibrate(int n_iter, int n_thin);
    DensityData get_predictive_density(
            int n_posterior_samples, int n_points, double quantile_edge_width);
    DensityOutput get_posterior_calendar_age_density(int output_offset, int ident);
    DensityOutput get_single_calendar_age_likelihood(int output_offset, int ident);

    std::vector<double> get_c14_age() { return c14_age; }
    std::vector<double> get_c14_sig() { return c14_sig; }

    double get_lambda()  { return lambda; }
    double get_nu1() { return nu1; }
    double get_nu2() { return nu2; }
    double get_A() { return A; }
    double get_B() { return B; }
    double get_alpha_shape() { return alpha_shape; }
    double get_alpha_rate() { return alpha_rate; }
    double get_slice_width() { return slice_width; }
    double get_slice_multiplier() { return slice_multiplier; }

    std::vector<double> get_alpha() { return alpha; }
    std::vector<int> get_n_clust() { return n_clust; }
    std::vector<double> get_mu_phi() { return mu_phi; }
    std::vector<double> get_phi() { return phi_i; }
    std::vector<double> get_tau() { return tau_i; }
    std::vector<double> get_weight() { return weight_i; }
    std::vector<double> get_calendar_age(int ident) {
        std::vector<double> output(n_out);
        for (int i = 0; i < n_out; i++) output[i] = calendar_age[i][ident];
        return output;
    }
    std::vector<int> get_cluster_ids() { return cluster_ids_i; }
};

#endif //CARBONDATE_WALKERDPMM_H
