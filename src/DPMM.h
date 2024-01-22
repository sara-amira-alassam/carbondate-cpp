/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_DPMM_H
#define CARBONDATE_DPMM_H
#include "carbondate_internal.h"
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

/*
 * This is a base class for one of to available parent classes: WalkerDPMM and PolyaUrnDPMM.
 * It contains all the public methods for these to classes needed to perform a calculation, however some of the
 * necessary private methods to implement these are be defined in the child classes (i.e. this class cannot be used
 * as-is to perform a complete calculation, instead you need an instance of WalkerDPMM or PolyaUrnDPMM).
 *
 * Note that it expects the global variable `project_name` to be present (which is achieved by call the read_arguments()
 * function.
 */
class DPMM {
protected:
    int n_work_update = 5000; // How many iterations between updating the work file

    std::vector<double> rc_determinations;  // observed radiocarbon determinations
    std::vector<double> rc_sigmas;  // radiocarbon determination uncertainties
    bool f14c_inputs{};     // Whether the radiocarbon determinations are c14 age or f14c age
    CalCurve calcurve;            // original calibration curve data
    YearlyCalCurve yearwise_calcurve;   // calibration curve interpolated for every year of calendar age

    int n_obs{}, n_out{};

    // Hyperparameters
    double lambda{}, nu1{}, nu2{};
    double A{}, B{};
    double alpha_shape = 1., alpha_rate = 1.;
    double slice_width{}, slice_multiplier = 10.;

    // Instant value of calendar age and stored values
    std::vector<double> calendar_age_i;
    std::vector<std::vector<double>> calendar_age;

    // Instant values of DPMM parameters
    // n_clust refers to the number of unique cluster ids for the observations
    int n_clust_i{};
    double alpha_i{}, mu_phi_i{};
    std::vector<double> phi_i, tau_i;
    std::vector<int> cluster_ids_i;

    // Stored values of DPMM parameters and output
    std::vector<int> n_clust;
    std::vector<double> alpha, mu_phi;
    std::vector<std::vector<double>> phi, tau;
    std::vector<std::vector<int>> cluster_ids;

    // SPD ranges - used for initialising certain start values and hyperparameters
    std::vector<double> spd_range_1_sigma = {0., 0.}, spd_range_2_sigma = {0., 0.}, spd_range_3_sigma = {0., 0.};

protected:
    virtual void _resize_storage();
    void _initialise_calendar_age_and_spd_ranges();
    void _initialise_hyperparameters();
    virtual void _initialise_clusters();
    void _interpolate_calibration_curve();
    virtual void _perform_update_step() {};
    virtual void _store_current_values(int i);
    double _cal_age_log_likelihood(
            double cal_age,
            double prmean,
            double prsig,
            double obs_c14_age,
            double obs_c14_sig);
    void _update_cluster_phi_and_tau(int cluster_id, const std::vector<double>& cluster_calendar_ages);
    void _update_mu_phi();
    void _update_alpha();
    double _alpha_log_prior(double alpha_value);
    virtual double _alpha_log_likelihood(double alpha_value);
    void _update_calendar_ages();
    double _log_marginal_normal_gamma(double cal_age, double mu_phi_s);
    virtual double _calculate_density_sample(int sample_id, double calendar_age_BP);

public:
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
