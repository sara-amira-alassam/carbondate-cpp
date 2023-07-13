#ifndef CARBONDATE_DENSITYOUTPUT_H
#define CARBONDATE_DENSITYOUTPUT_H
#include "carbondate.h"

class UnableToWriteToOutputFileException : public CarbondateException {
public:
    explicit UnableToWriteToOutputFileException(const std::string& file_path) {
        _error_message = "Unable to write to output file " + file_path;
    }
};


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
    std::string variable_line(const std::string& var_name, const std::string& var);
    std::string output_line(const std::string& var_name, double var);
    std::string output_line(const std::string& var_name, const std::vector<double>& var);
    std::string comment_line(const std::string &comment, int &comment_index);

protected:
    virtual std::vector<std::string> get_output_lines();
    virtual std::string range_lines(int &comment_index);

public:
    DensityOutput(int index, double resolution);
    void set_probability(const std::vector<double>& probability);
    void print();
};


#endif //CARBONDATE_DENSITYOUTPUT_H
