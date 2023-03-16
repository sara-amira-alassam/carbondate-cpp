//
// Created by Sara Al-Assam on 21/02/2023.
//

#ifndef CARBONDATE_DENSITYOUTPUT_H
#define CARBONDATE_DENSITYOUTPUT_H
#include <vector>
#include <string>

class DensityOutput {
    std::string _output_var;
    std::string _output_prefix;
    std::vector<double> _prob_yearwise;
    double _prob_max = 0.;
    double _date;
    double _error;
    std::string _name;

    std::vector<double> _prob_smoothed;
    int _resolution_smoothed{};
    double _start_calAD_smoothed{};
    double _prob_norm_smoothed{};

public:
    double start_calAD = 0;
    double mean_calAD = 0;
    double sigma = 0;
    double median_calAD = 0;

private:
    std::string variable_line(const std::string& var_name, int var);
    std::string variable_line(const std::string& var_name, double var);
    std::string variable_line(const std::string& var_name, const std::string& var);
    std::string output_line(const std::string& var_name, int var);
    std::string output_line(const std::string& var_name, double var);
    std::string output_line(const std::string& var_name, const std::vector<double>& var);
    std::string comment_line(const std::string &comment, int &comment_index);
    std::string range_lines(
            int range_index, double probability, int resolution, int& comment_index);

    void calculate_probability_smoothed(int resolution);
    double find_probability_and_ranges_for_cut_off_smoothed(
            double cut_off, std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> get_ranges_by_bisection(double probability, int resolution);
    std::vector<std::string> get_output_lines(int resolution);

    // Get rid of the below functions after testing
/*    double find_probability_and_ranges_for_cut_off(
            double cut_off, std::vector<std::vector<double>>& ranges);
    std::vector<std::vector<double>> get_ranges_by_bisection(double probability);
    std::vector<std::vector<double>> get_ranges(double probability);*/

public:
    DensityOutput(
            double date,
            double error,
            const std::string& name);
    void set_yearwise_probability(std::vector<double> probability);
    void write_to_file(
            int resolution,
            const std::string& file_prefix,
            const std::string& output_var,
            const std::string& output_name);
    std::string get_name() { return _name; }
    double get_lower() { return start_calAD; }
    double get_upper() { return start_calAD + (double) _prob_yearwise.size(); }
};


#endif //CARBONDATE_DENSITYOUTPUT_H
