#ifndef CARBONDATE_PLAIN_TEXT_H
#define CARBONDATE_PLAIN_TEXT_H

#include "carbondate.h"

class UnableToFindTextFileException : public CarbondateException {
public:
    explicit UnableToFindTextFileException(const std::string& file_path) {
        _error_message = "Unable to find text file at " + file_path;
    }
};

class UnableToWriteToTextFileException : public CarbondateException {
public:
    explicit UnableToWriteToTextFileException(const std::string& file_path) {
        _error_message = "Unable to write to text file at " + file_path;
    }
};

void initialize_text_file();
void update_text_file(const std::string& label, const std::vector<double>& values);

#endif //CARBONDATE_PLAIN_TEXT_H
