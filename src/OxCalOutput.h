#ifndef CARBONDATE_OXCALOUTPUT_H
#define CARBONDATE_OXCALOUTPUT_H

#include "DensityOutput.h"

class OxCalOutput {
    std::vector<DensityOutput> _likelihoods;
    std::vector<DensityOutput> _posteriors;
    int _resolution;
    std::string _file_prefix;

private:
    void initialise_file();
    void print_calibration_data();
    void print_model();
    void print_likelihoods();
    void print_predictive_density();
    void print_posteriors();

public:
    OxCalOutput(int n_obs, int resolution, const std::string& file_prefix);
    void set_likelihood(int ident, const DensityOutput& likelihood);
    void set_posterior(int ident, const DensityOutput& posterior);
    void print();

};


#endif //CARBONDATE_OXCALOUTPUT_H
