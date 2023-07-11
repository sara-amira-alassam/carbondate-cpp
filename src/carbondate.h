#ifndef CARBONDATE_CARBONDATE_H
#define CARBONDATE_CARBONDATE_H

#ifndef CARBONDATE_VERSION
#define CARBONDATE_VERSION "1.0.0"
#endif

#include <exception>
#include <string>
#include <vector>

#include "helpers.h"

class CarbondateException: public std::exception {
public:
    const char* what() const noexcept override {
        return _error_message.c_str();
    }
protected:
    std::string _error_message;
};

#endif //CARBONDATE_CARBONDATE_H
