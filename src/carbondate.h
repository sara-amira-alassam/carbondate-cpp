#ifndef CARBONDATE_CARBONDATE_INTERNAL_H
#define CARBONDATE_CARBONDATE_INTERNAL_H

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

inline std::string carbondate_short_reference() {
    return "Carbondate v" + (std::string) CARBONDATE_VERSION;
}

inline std::string carbondate_long_reference() {
    return "Carbondate v" + (std::string) CARBONDATE_VERSION + " TJ Heaton (2023)";
}
#endif //CARBONDATE_CARBONDATE_INTERNAL_H
