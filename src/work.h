#ifndef CARBONDATE_WORK_H
#define CARBONDATE_WORK_H

#include "carbondate.h"

class WorkFileRemovedException : public CarbondateException {
private:
    std::string _error_message = "Work file removed whilst program is running";
};

class UnableToCreateWorkFileException : public CarbondateException {
private:
    std::string _error_message = "Could not create a writeable work file";
};

class UnableToWriteToWorkFileException : public CarbondateException {
private:
    std::string _error_message = "Could not write to the work file";
};

void create_work_file();

void update_work_file_mcmc(double done, int iterations);

void update_work_file_postprocessing(int iterations);

void check_for_work_file();

void remove_work_file();

#endif //CARBONDATE_WORK_H
