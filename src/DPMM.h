//
// Created by Sara Admin on 21/06/2023.
//

#ifndef CARBONDATE_DPMM_H
#define CARBONDATE_DPMM_H
#include <utility>
#include <vector>
#include <string>
#define MATHLIB_STANDALONE
#include <Rmath.h>

struct CalCurve {
    std::vector<double> cal_age;
    std::vector<double> c14_age;
    std::vector<double> c14_sig;
    std::vector<double> f14c_age;
    std::vector<double> f14c_sig;
};

struct YearlyCalCurve {
    std::vector<double> cal_age;
    std::vector<double> rc_age;
    std::vector<double> rc_sig;
};

struct DensityData {
    std::vector<double> cal_age_AD;
    std::vector<double> mean;
    std::vector<double> ci_lower;
    std::vector<double> ci_upper;

    explicit DensityData(int n_points)
            : cal_age_AD(n_points), mean(n_points), ci_lower(n_points), ci_upper(n_points) {};
};

class DPMM {
protected:
    std::string _file_prefix;

    std::vector<double> rc_determinations;  // observed radiocarbon determinations
    std::vector<double> rc_sigmas;  // radiocarbon determination uncertainties
    bool f14c_inputs;     // Whether the radiocarbon determinations are c14 age or f14c age
    CalCurve calcurve;            // original calibration curve data
    YearlyCalCurve yearwise_calcurve;   // calibration curve interpolated for every year of calendar age

    int n_obs, n_out;

    // Hyperparameters
    double lambda, nu1, nu2;
    double A, B;
    double alpha_shape = 1., alpha_rate = 1.;
    double slice_width, slice_multiplier = 10.;

    // Instant value of calendar age and stored values
    std::vector<double> calendar_age_i;
    std::vector<std::vector<double>> calendar_age;

    // Instant values of DPMM parameters
    // n_clust refers to the number of unique cluster ids for the observations
    int n_clust_i;
    double alpha_i, mu_phi_i;
    std::vector<double> phi_i, tau_i;
    std::vector<int> cluster_ids_i;

    // Stored values of DPMM parameters and output
    std::vector<int> n_clust;
    std::vector<double> alpha, mu_phi;
    std::vector<std::vector<double>> phi, tau;
    std::vector<std::vector<int>> cluster_ids;

    // SPD ranges - used for initialising certain value
    std::vector<double> spd_range_1_sigma = {0., 0.}, spd_range_2_sigma = {0., 0.}, spd_range_3_sigma = {0., 0.};

protected:
    virtual void initialise_storage();
    void initialise_calendar_age_and_spd_ranges();
    void initialise_hyperparameters();
    virtual void initialise_clusters();
    void interpolate_calibration_curve();
    virtual void perform_update_step() {};
    virtual void store_current_values(int i);
    double cal_age_log_likelihood(
            double cal_age,
            double prmean,
            double prsig,
            double obs_c14_age,
            double obs_c14_sig);
    void update_cluster_phi_and_tau(int cluster_id, const std::vector<double>& cluster_calendar_ages);
    void update_mu_phi();
    void update_alpha();
    double alpha_log_prior(double alpha_value);
    virtual double alpha_log_likelihood(double alpha_value);
    void update_calendar_ages();
    double log_marginal_normal_gamma(double cal_age, double mu_phi_s);
    virtual double calculate_density_sample(int sample_id, double calendar_age_BP);

public:
    DPMM(std::string file_prefix): _file_prefix(std::move(file_prefix)) {}
    void initialise(
            std::vector<double> i_rc_determinations,
            std::vector<double> i_rc_sigmas,
            bool i_f14c_inputs,
            std::vector<double> cc_cal_age,
            std::vector<double> cc_c14_age,
            std::vector<double> cc_c14_sig,
            int rng_seed = 0);
    void calibrate(int n_iter, int n_thin);
    DensityData get_predictive_density(int n_posterior_samples, double resolution, double quantile_edge_width);
    std::vector<double> get_posterior_calendar_ages(int ident);
    int get_nobs() const { return n_obs; }
};


#endif //CARBONDATE_DPMM_H
