#include "work.h"
#include <cstdio>
#include <fstream>

void create_work_file(const std::string& file_prefix){
    std::string filepath = file_prefix + ".work";

    std::ofstream file(filepath);
    if (file.is_open()) {
        file.close();
    } else {
        throw UnableToCreateWorkFileException();
    }
}

void update_work_file_mcmc(const std::string& file_prefix, double done, int passes){

}

void update_work_file_postprocessing(const std::string& file_prefix){

}

void check_for_work_file(const std::string& file_prefix){
    std::string filepath = file_prefix + ".work";
    std::ifstream file(filepath.c_str());
    if (!file.good()) {
        throw WorkFileRemovedException();
    }
}

void remove_work_file(const std::string& file_prefix) {
    std::string filepath = file_prefix + ".work";

    std::remove(filepath.c_str());
}
