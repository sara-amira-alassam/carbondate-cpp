//
// Created by Sara Al-Assam on 20/03/2023.
//

#ifndef CARBONDATE_POSTERIORDENSITYOUTPUT_H
#define CARBONDATE_POSTERIORDENSITYOUTPUT_H

#include "DensityOutput.h"

class PosteriorDensityOutput : public DensityOutput {
public:
    PosteriorDensityOutput(
            int ident,
            int offset,
            double resolution,
            const std::vector<bool>& ranges,
            const std::vector<double>& posterior_calendar_ages_AD);
};


#endif //CARBONDATE_POSTERIORDENSITYOUTPUT_H
