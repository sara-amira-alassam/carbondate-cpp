//
// Created by Sara Al-Assam on 20/03/2023.
//

#include "PosteriorDensityOutput.h"

PosteriorDensityOutput::PosteriorDensityOutput(int ident, double resolution)
        : DensityOutput(ident + 2, resolution) {
    _output_var = "ocd[" + std::to_string(_index + 1) + "]";
    _output_prefix = _output_var + ".posterior";
}
