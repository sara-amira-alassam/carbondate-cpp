#include <algorithm>
#include <string>
/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "helpers.h"

// Mean of elements in a vector
double mean(const std::vector<double>& vec) {
    double mean = 0.0;
    for (double elem : vec) mean += elem;
    mean /= (double) vec.size();
    return mean;
}

// Mean of a probability distribution, assumes normalised
double mean(const std::vector<double>& vec, const std::vector<double>& probability) {
    double mean = 0.;
    for (int i = 0; i < vec.size(); i++) mean += vec[i] * probability[i];
    return mean;
}

// Standard deviation of elements in a vector
double sigma(const std::vector<double>& vec, double mean) {
    double sigma_squared = 0.;
    for (double elem : vec) sigma_squared += pow(elem - mean, 2);
    return sqrt(sigma_squared / (double) vec.size());
}

// Standard deviation of a probability distribution, assumes normalised
double sigma(const std::vector<double>& vec, const std::vector<double>& probability, double mean) {
    double sigma_squared = 0.;
    for (int i = 0; i < vec.size(); i++) sigma_squared += pow(vec[i] - mean, 2) * probability[i];
    return sqrt(sigma_squared);
}

// Middle value of sorted elements in a vector
double median(std::vector <double> vec) {
   unsigned n = vec.size(), half = n/2;

    if (n % 2 == 0) {
        std::nth_element(vec.begin(), vec.begin() + half - 1, vec.end());
        std::nth_element(vec.begin() + half, vec.begin() + half, vec.end());
        return (vec[half - 1] + vec[half])/2.;
    }
    std::nth_element(vec.begin(), vec.begin() + half, vec.end());
    return vec[half];
}

// Median of a probability distribution, assumes normalised
double median(const std::vector<double>& vec, const std::vector<double>& probability) {
    double cumulative_probability = 0;
    for (int i = 0; i < vec.size(); i++) {
        cumulative_probability += probability[i];
        if (cumulative_probability >= 0.5) return vec[i];
    }
    return 0.;
}

double mad(std::vector <double> vec) {
    double constant = 1.4826, centre = median(vec);
    unsigned n = vec.size();
    std::vector <double> deviation(n);

    for (int i = 0; i < n; i++) deviation[i] = abs(vec[i] - centre);

    return constant * median(deviation);
}

// Returns the value of the maximum element in a vector - the minimum element in a vector
double max_diff(std::vector<double> vec) {
    auto minmax = std::minmax_element(vec.begin(), vec.end());
    return *(minmax.second) - *(minmax.first);
}

// Finds the point yi that corresponds to xi, given known points (x1, y1) and (x2, y2)
double interpolate_linear(double xi, double x1, double x2, double y1, double y2) {
    return y1 + (y2 - y1) * ((xi - x1)/(x2 - x1));
}

// Finds the index for which all values below this are less than the cut-off
int get_left_boundary(const std::vector<double>& vec, double cut_off) {
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] >= cut_off) return i;
    }
    return 0;
}

// Finds the index for which all values above this are more than the cut-off
int get_right_boundary(const std::vector<double>& vec, double cut_off) {
    for (int i = (int) vec.size() - 1; i >= 0; i--) {
        if (vec[i] <= cut_off) return i;
    }
    return (int) vec.size() - 1;
}

// Finds the quantiles at a given edge width away from the start and end of the distribution,
// where the distribution is defined by a vector of sampled values
// Note this partially sorts the vector provided as an argument
void edge_quantiles(
        std::vector<double>& vec,   // Vector to find quantiles of. Assumed unsorted
        double edge_width,         // Probability distance from each side of the distribution
        double& lower_quantile,    // Lower quantile value to return
        double& upper_quantile) {  // Upper quantile value to return

    double indl, indu, gl, gu;
    int jl, ju;

    // The indices we use are defined by type = 7 for the R quantile function
    indl = edge_width * ((double) vec.size() - 1.) + 1.;
    jl = std::floor(indl);
    indu = (1. - edge_width) * ((double) vec.size() - 1.) + 1.;
    ju = std::floor(indu);

    // Rather than sorting the entire vector, just sort partially to find the elements we'll look up
    std::nth_element(vec.begin(), vec.begin() + jl - 1, vec.end());
    std::nth_element(vec.begin() + jl, vec.begin() + jl, vec.end());
    std::nth_element(vec.begin() + jl + 1, vec.begin() + ju - 1, vec.end());
    std::nth_element(vec.begin() + ju, vec.begin() + ju, vec.end());

    // quantiles found using the formula for type = 7 in the R quantile function
    gl = indl - jl;
    lower_quantile = (1. - gl) * vec[jl - 1] + gl * vec[jl];

    gu = indu - ju;
    upper_quantile = (1. - gu) * vec[ju - 1] + gu * vec[ju];
}

// Adapted from do_sample in R/src/main/sample.c and from EmpiricalSample in Rcpp package
void get_sample_ids(std::vector<int>& ans, int start_index, int finish_index) {

    int n = finish_index - start_index + 1;
    bool replace = ans.size() >= n;

    if (replace || ans.size() < 2) {
        for (int & an : ans) an = static_cast<int>(R_unif_index(n)) + start_index;
        return;
    }

    std::vector<int> x(n);
    for (int i = 0; i < n; i++) x[i] = i;

    for (int & an : ans) {
        int j = static_cast<int>(R_unif_index(n));
        an = x[j] + start_index;
        x[j] = x[--n];
    }
}

// Adapted from `ProbSampleNoReplace` in R/src/main/sample.c and from SampleNoReplace
// in Rcpp package, but substantially simplified since we know we're only ever
// going to be calling this with sz = 1. Tested that it gives the same result as
// calling sample.int from R.
int sample_integer(unsigned n, std::vector<double> prob, bool one_based) {

    int adj = one_based ? 0 : 1;
    double rT, mass, sum_p = 0.;
    int i, j;
    std::vector<double> p(n);
    std::vector<int> perm(n);

    for (i = 0; i < n; i++) {
        perm[i] = i + 1;
        if (R_FINITE(prob[i]) && prob[i] > 0.0) {
            p[i] = prob[i];
            sum_p += p[i];
        } else {
            p[i] = 0.0;
        }
    }
    revsort(&p[0], &perm[0], (int) n);

    rT = unif_rand() * sum_p;
    mass = 0.0;
    for (j = 0; j < n-1; j++) {
        mass += p[j];
        if (rT <= mass) break;
    }
    return perm[j] - adj;
}

void update_progress_bar(double progress) {
    int bar_width = 100;
    std::string bar_string(bar_width, '=');

    int val = (int) (progress * 100);
    int left_padding = (int) (progress * bar_width);
    int right_padding = bar_width - left_padding;
    printf("\r%3d%% [%.*s%*s]", val, left_padding, bar_string.c_str(), right_padding, "");
    fflush(stdout);
}

double to_calAD(double year_calPB) {
    return 1950.5 - year_calPB;
}

std::string to_string(double var, int max_digits) {
    std::string string_var, format_string = "%." + std::to_string(max_digits) + "g";
    char temp_string[32] = "";
    snprintf(temp_string, 13, format_string.c_str(), var);
    string_var = temp_string;
    return string_var;
}

std::string to_percent_string(double fraction) {
    std::string percent;
    char temp_string[8] = "";
    snprintf(temp_string, 8, "%4.1f%%", fraction * 100.);
    percent = temp_string;
    return percent;
}

void convert_to_c14_age(
        const std::vector<double>&f14c_age,
        const std::vector<double>&f14c_sig,
        std::vector<double>&c14_age,
        std::vector<double>&c14_sig) {

    unsigned n = f14c_age.size();

    c14_age.resize(n);
    c14_sig.resize(n);
    for (int i = 0; i < n; i++) {
        c14_age[i] = -8033. * log(f14c_age[i]);
        c14_sig[i] =  8033. * f14c_sig[i] / f14c_age[i];
    }
}

void convert_to_f14c_age(
        const std::vector<double>&c14_age,
        const std::vector<double>&c14_sig,
        std::vector<double>&f14c_age,
        std::vector<double>&f14c_sig){

    unsigned n = c14_age.size();

    f14c_age.resize(n);
    f14c_sig.resize(n);
    for (int i = 0; i < n; i++) {
        f14c_age[i] = exp(c14_age[i] / -8033.);
        f14c_sig[i] = f14c_age[i] * c14_sig[i] / 8033.;
    }
}
