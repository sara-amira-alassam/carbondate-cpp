#ifndef CARBONDATE_LOG_H
#define CARBONDATE_LOG_H

#include "carbondate.h"


class UnableToWriteToLogFileException : public CarbondateException {
public:
    explicit UnableToWriteToLogFileException(const std::string& file_path) {
        _error_message = "Unable to write to log file at " + file_path;
    }
};

void initialize_log_file();
void update_log_file(const std::string& log_line);

#endif //CARBONDATE_LOG_H
