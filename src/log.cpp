#include <fstream>
#include "log.h"

std::string log_file_path() {
    return project_directory + project_name + ".log";
}

/*
 * Must be called before any calls to update the log file. Checks it exists and adds the carbondate reference.
 * Note there must first be a call to read_arguments() to set the prefix for the file.
 */
void initialize_log_file() {
    std::string filepath = log_file_path();

    update_log_file(carbondate_long_reference());
}

/* Adds line to the log file. Note initialize_log_file() must be called once before any call to update the log file */
void update_log_file(const std::string& log_line) {
    std::string filepath = log_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToLogFileException(filepath);

    file << log_line << std::endl;
    file.close();
}
