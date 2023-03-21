//
// Created by Sara Al-Assam on 21/02/2023.
//

#ifndef CARBONDATE_DENSITYOUTPUT_H
#define CARBONDATE_DENSITYOUTPUT_H
#include <vector>
#include <string>

class DensityOutput {
    std::vector<double> _probability;
    double _prob_norm = 0.;

protected:
    std::string _output_var;
    std::string _output_prefix;
    int _index;
    double _resolution;
    double _prob_max = 0.;

public:
    double start_calAD = 0;
    double mean_calAD = 0;
    double sigma = 0;
    double median_calAD = 0;

private:
    std::string variable_line(const std::string& var_name, const std::string& var);
    std::string output_line(const std::string& var_name, int var);
    std::string output_line(const std::string& var_name, double var);
    std::string output_line(const std::string& var_name, const std::vector<double>& var);
    std::string comment_line(const std::string &comment, int &comment_index);
    std::string range_lines(int range_index, double probability, int& comment_index);

    double find_probability_and_ranges_for_cut_off(
            double cut_off, std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> get_ranges_by_bisection(double probability);

protected:
    virtual std::vector<std::string> get_output_lines();

public:
    DensityOutput(int index, double resolution);

    void set_probability(const std::vector<double>& probability);
    void write_to_file(
            const std::string& file_prefix,
            const std::string& output_var,
            const std::string& output_name);
};


#endif //CARBONDATE_DENSITYOUTPUT_H
