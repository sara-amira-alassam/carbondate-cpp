#ifndef CARBONDATE_OXCALOUTPUT_H
#define CARBONDATE_OXCALOUTPUT_H

#include "DensityOutput.h"

class OxCalOutput {
    std::vector<DensityOutput> _posteriors;
    std::string _file_prefix;

private:
    void initialise_file();
    void print_model();
    void print_predictive_density();
    void print_posteriors();
    void append_to_file(const std::vector<std::string>& output_lines);

public:
    DensityOutput predictive_density;
    OxCalOutput(int n_obs, const std::string& file_prefix);
    void set_posterior(int ident, const DensityOutput& posterior);
    void print();

};


#endif //CARBONDATE_OXCALOUTPUT_H
