#ifndef CARBONDATE_WORK_H
#define CARBONDATE_WORK_H

#include "carbondate_internal.h"

class WorkFileRemovedException : public CarbondateException {
public:
    explicit WorkFileRemovedException(const std::string& file_path) {
        _error_message = "Work file at " + file_path + " removed whilst program is running";
    }
};

class UnableToCreateWorkFileException : public CarbondateException {
public:
    explicit UnableToCreateWorkFileException(const std::string& file_path) {
        _error_message = "Unable to create writable work file at " + file_path;
    }
};

class UnableToWriteToWorkFileException : public CarbondateException {
public:
    explicit UnableToWriteToWorkFileException(const std::string& file_path) {
        _error_message = "Unable to write to work file at " + file_path;
    }
};

void create_work_file();

void update_work_file_mcmc(double done, int iterations);

void update_work_file_postprocessing(int iterations);

void check_for_work_file();

void remove_work_file();

#endif //CARBONDATE_WORK_H
