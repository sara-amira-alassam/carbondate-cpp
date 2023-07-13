#ifndef CARBONDATE_LOG_H
#define CARBONDATE_LOG_H

#include "carbondate.h"

class UnableToFindLogFileException : public CarbondateException {
public:
    explicit UnableToFindLogFileException(const std::string& file_path) {
        _error_message = "Unable to find log file at " + file_path;
    }
};

class UnableToWriteToLogFileException : public CarbondateException {
public:
    explicit UnableToWriteToLogFileException(const std::string& file_path) {
        _error_message = "Unable to write to log file at " + file_path;
    }
};

void initialize_log_file(const std::string& file_prefix);
void update_log_file(const std::string& file_prefix, const std::string& log_line);

#endif //CARBONDATE_LOG_H
