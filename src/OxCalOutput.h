#ifndef CARBONDATE_OXCALOUTPUT_H
#define CARBONDATE_OXCALOUTPUT_H

#include "PosteriorDensityOutput.h"
#include "PredictiveDensityOutput.h"

class OxCalOutput {
    std::vector<PosteriorDensityOutput> _posteriors;
    std::string _file_prefix;

private:
    void initialise_file();
    void print_model();
    void print_predictive_density();
    void print_posteriors();
    void append_to_file(const std::vector<std::string>& output_lines);

public:
    PredictiveDensityOutput _predictive_density;
    OxCalOutput(
            int n_obs, const std::string &file_prefix, PredictiveDensityOutput predictiveDensity);
    void set_posterior(int ident, const PosteriorDensityOutput& posterior);
    void print();

};


#endif //CARBONDATE_OXCALOUTPUT_H
