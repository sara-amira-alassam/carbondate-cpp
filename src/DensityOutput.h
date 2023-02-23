//
// Created by Sara Al-Assam on 21/02/2023.
//

#ifndef CARBONDATE_DENSITYOUTPUT_H
#define CARBONDATE_DENSITYOUTPUT_H
#include <vector>
#include <string>

class DensityOutput {
    std::string _output_var;
    int _index;
    std::string _output_name;
    std::string _output_prefix;

public:
    std::vector<double> prob;
    double prob_norm = 0;
    double start_calAD = 0;
    int resolution = 0;
    double mean_calAD = 0;
    double sigma = 0;
    double median_calAD = 0;

private:
    std::string output_prefix();
    std::string output_line(const std::string& var_name, int var);
    std::string output_line(const std::string& var_name, double var);
    std::string output_line(const std::string& var_name, const std::vector<double>& var);

public:
    DensityOutput(std::string output_var, int index, std::string output_name);
    void print();
    std::vector<std::vector<double>> as_columns();
};


#endif //CARBONDATE_DENSITYOUTPUT_H
