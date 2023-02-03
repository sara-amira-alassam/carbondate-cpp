#include <vector>
#include <algorithm>
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "helpers.h"


double mean(const std::vector<double>& vec) {
    double mean = 0.0;
    for (double elem : vec) {
        mean += elem;
    }
    mean /= (double) vec.size();
    return mean;
}

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

double mad(std::vector <double> vec) {
    double constant = 1.4826, centre = median(vec);
    unsigned n = vec.size();
    std::vector <double> deviation(n);

    for (int i = 0; i < n; i++) {
        deviation[i] = abs(vec[i] - centre);
    }
    return constant * median(deviation);
}

double max_diff(std::vector<double> vec) {
    auto minmax = std::minmax_element(vec.begin(), vec.end());
    return *(minmax.second) - *(minmax.first);
}

double interpolate_linear(double xi, double x1, double x2, double y1, double y2) {
    return y1 + (y2 - y1) * ((xi - x1)/(x2 - x1));
}

// Finds the quantiles at a given edge width away from the start and end of the distribution.
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
void get_sample_ids(std::vector<int>& ans, int start_index, int finish_index,int size) {

    ans.resize(size);
    int n = finish_index - start_index + 1;
    bool replace = size >= n;

    if (replace || size < 2) {
        for (int i = 0 ; i < size; i++) {
            ans[i] = static_cast<int>(R_unif_index(n)) + start_index;
        }
        return;
    }

    std::vector<int> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = i;
    }

    for (int i = 0 ; i < size; i++) {
        int j = static_cast<int>(R_unif_index(n));
        ans[i] = x[j] + start_index;
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
