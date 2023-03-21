#ifndef CARBONDATE_OXCALOUTPUT_H
#define CARBONDATE_OXCALOUTPUT_H

#include "PosteriorDensityOutput.h"
#include "PredictiveDensityOutput.h"

class OxCalOutput {
    std::string _file_prefix;

private:
    void initialise_file();

public:
    OxCalOutput(std::string file_prefix);
};

#endif //CARBONDATE_OXCALOUTPUT_H
