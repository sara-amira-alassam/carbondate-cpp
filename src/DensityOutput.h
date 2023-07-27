#ifndef CARBONDATE_DENSITYOUTPUT_H
#define CARBONDATE_DENSITYOUTPUT_H
#include "carbondate.h"

class UnableToWriteToOutputFileException : public CarbondateException {
public:
    explicit UnableToWriteToOutputFileException(const std::string& file_path) {
        _error_message = "Unable to write to output file " + file_path;
    }
};

class DensityOutputException : public CarbondateException {
public:
    explicit DensityOutputException(const std::string& error_message) {
        _error_message = error_message;
    }
};

/*
 * This class is used to create density output in a format that is recognised by the OxCal front-end, using the only
 * public method print() which writes the the *.js file.
 *
 * Note that it expects the global variable `project_name` to be present (which is achieved by call the read_arguments()
 * function.
 *
 * Note this is a base class for the two supported types of output - the predictive density and the individual posterior
 * densities - and the child classes should be used to create output as there are specific overridden methods on each
 * for writing the output lines.
 */
class DensityOutput {
protected:
    std::vector<double> _probability;
    std::string _output_var;
    std::string _output_prefix;
    int _index;
    double _resolution;
    double _prob_max = 0.;
    double _prob_norm = 0.;
    double _start_calAD = 0;
    double _mean_calAD = 0;
    double _sigma_calAD = 0;
    double _median_calAD = 0;

protected:
    std::string _variable_line(const std::string& var_name, const std::string& var);
    std::string _output_line(const std::string& var_name, double var);
    std::string _output_line(const std::string& var_name, const std::vector<double>& var);
    std::string _comment_line(const std::string &comment, int &comment_index);
    void _set_probability(const std::vector<double>& probability);

protected:
    virtual std::vector<std::string> _get_output_lines();
    virtual std::string _range_lines(int &comment_index);

public:
    DensityOutput(int index, double resolution);
    void print();
};


#endif //CARBONDATE_DENSITYOUTPUT_H
