#include <fstream>
#include "log.h"

#ifndef LOG_PREFIX
#define LOG_PREFIX "../output/"
#endif

std::string log_file_path(const std::string& file_prefix) {
    return LOG_PREFIX + file_prefix + ".log";
}

void initialize_log_file(const std::string& file_prefix) {
    std::string filepath = log_file_path(file_prefix);

    std::ofstream file(filepath, std::fstream::in);
    if (!file.is_open()) {
        throw UnableToFindLogFileException(filepath);
    }
    file.close();
    update_log_file(file_prefix, carbondate_long_reference());
}

void update_log_file(const std::string& file_prefix, const std::string& log_line) {
    std::string filepath = log_file_path(file_prefix);

    std::ofstream file(filepath, std::fstream::app);
    if (file.is_open()) {
        file << log_line << std::endl;
        file.close();
    } else {
        throw UnableToWriteToLogFileException(filepath);
    }
}
