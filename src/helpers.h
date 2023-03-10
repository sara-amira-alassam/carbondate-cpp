//
// Created by Sara Al-Assam on 02/02/2023.
//
#include <vector>
#ifndef CARBONDATE_HELPERS_H
#define CARBONDATE_HELPERS_H

double mean(const std::vector<double>& vec);
double mean(const std::vector<double>& vec, const std::vector<double>& probability);
double sigma(const std::vector<double>& vec, double mean);
double sigma(const std::vector<double>& vec, const std::vector<double>& probability, double mean);
double median(std::vector <double> vec);
double median(const std::vector<double>& vec, const std::vector<double>& probability);
double mad(std::vector <double> vec);
double max_diff(std::vector <double> vec);
double interpolate_linear(double xi, double x1, double x2, double y1, double y2);
int get_left_boundary(const std::vector<double>& vec, double cut_off);
int get_right_boundary(const std::vector<double>& vec, double cut_off);
void edge_quantiles(
        std::vector<double>& vec, double edge_width, double& lower_quantile, double& upper_quantile);
void get_sample_ids(std::vector<int>& ans, int start_index, int finish_index,int size);
int sample_integer(unsigned n, std::vector<double> prob, bool one_based);
void rsort_with_index(double *x, int *indx,int n);
void revsort(double *a, int *ib, int n);
void update_progress_bar(double progress);
double to_calAD(double year_calPB);

#endif //CARBONDATE_HELPERS_H
