#ifndef CARBONDATE_WORK_H
#define CARBONDATE_WORK_H

#include "carbondate.h"

class WorkFileRemovedException : public CarbondateException {
public:
    WorkFileRemovedException() {
        _error_message = "Work file removed whilst program is running";
    }
};

class UnableToCreateWorkFileException : public CarbondateException {
private:
    std::string _error_message = "Could not create a writeable work file";
};

class UnableToWriteToWorkFileException : public CarbondateException {
private:
    std::string _error_message = "Could not write to the work file";
};

void create_work_file(const std::string& file_prefix);

void update_work_file_mcmc(const std::string& file_prefix, double done, int iterations);

void update_work_file_postprocessing(const std::string& file_prefix, int iterations);

void check_for_work_file(const std::string& file_prefix);

void remove_work_file(const std::string& file_prefix);

#endif //CARBONDATE_WORK_H
