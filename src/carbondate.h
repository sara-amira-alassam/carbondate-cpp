#ifndef CARBONDATE_CARBONDATE_H
#define CARBONDATE_CARBONDATE_H

#include <exception>
#include <string>
#include <vector>

class CarbondateException: public std::exception {
public:
    const char* what() const noexcept override {
        return _error_message.c_str();
    }
protected:
    std::string _error_message;
};

#endif //CARBONDATE_CARBONDATE_H
