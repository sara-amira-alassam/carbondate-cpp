#ifndef CARBONDATE_CARBONDATE_INTERNAL_H
#define CARBONDATE_CARBONDATE_INTERNAL_H

#ifndef CARBONDATE_VERSION
#define CARBONDATE_VERSION "1.0.0"
#endif

#include <exception>
#include <string>
#include "helpers.h"

extern std::string project_name, project_directory, oxcal_version;

class CarbondateException: public std::exception {
public:
    const char* what() const noexcept override {
        return _error_message.c_str();
    }
protected:
    std::string _error_message;
};

inline std::string carbondate_short_reference() {
    return "carbondate v" + (std::string) CARBONDATE_VERSION;
}

inline std::string carbondate_full_reference() {
    return "carbondate v" + (std::string) CARBONDATE_VERSION + " Heaton et al (2024)";
}

inline std::string carbondate_long_reference() {
    return oxcal_version + "; carbondate v" + (std::string) CARBONDATE_VERSION + " Heaton (2024)";
}
#endif //CARBONDATE_CARBONDATE_INTERNAL_H
