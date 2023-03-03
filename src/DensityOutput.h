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
    std::vector<double> _prob_yearwise;
    double _prob_total = 0.;
    double _prob_max = 0.;

    std::vector<double> _prob_smoothed;
    int _resolution_smoothed;
    double _start_calAD_smoothed;
    double _prob_norm_smoothed;

public:
    double start_calAD = 0;
    double mean_calAD = 0;
    double sigma = 0;
    double median_calAD = 0;

private:
    std::string output_prefix();
    std::string output_line(const std::string& var_name, int var);
    std::string output_line(const std::string& var_name, double var);
    std::string output_line(const std::string& var_name, const std::vector<double>& var);
    void calculate_probability_smoothed(int resolution);
    double find_probability_and_ranges_for_cut_off(
            double cut_off, std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> get_ranges(double probability, int resolution);

public:
    DensityOutput(std::string output_var, int index, std::string output_name);
    void set_yearwise_probability(std::vector<double> probability);
    void print(int resolution);
    std::vector<std::vector<double>> as_columns(int resolution);
};


#endif //CARBONDATE_DENSITYOUTPUT_H
