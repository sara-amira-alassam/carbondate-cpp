#include <fstream>
#include "log.h"

#ifndef LOG_PREFIX
#define LOG_PREFIX "../output/"
#endif

std::string log_file_path() {
    return LOG_PREFIX + project_name + ".log";
}

void initialize_log_file() {
    std::string filepath = log_file_path();

    std::ofstream file(filepath, std::fstream::in);
    if (!file.is_open()) {
        throw UnableToFindLogFileException(filepath);
    }
    file.close();
    update_log_file(carbondate_long_reference());
}

void update_log_file(const std::string& log_line) {
    std::string filepath = log_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (file.is_open()) {
        file << log_line << std::endl;
        file.close();
    } else {
        throw UnableToWriteToLogFileException(filepath);
    }
}
