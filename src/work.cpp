#include "work.h"
#include <cstdio>
#include <fstream>

void create_work_file(const std::string& file_prefix) {
    std::string filepath = file_prefix + ".work";

    std::ofstream file(filepath);
    if (file.is_open()) {
        file.close();
    } else {
        throw UnableToCreateWorkFileException();
    }
}

void update_work_file_mcmc(const std::string& file_prefix, double done, int iterations) {
    std::string filepath = file_prefix + ".work";

    std::ofstream file(filepath, std::ios::trunc);
    if (file.is_open()) {
        std::string work_lines = "work.program=\"" + carbondate_short_reference() + "\";\n";
        work_lines += "work.operation=\"MCMC\";\n";
        work_lines += "work.done=" + to_string(done * 100., 3) + "; work.passes=" + std::to_string(iterations) + "; ";
        work_lines += "work.ok=100; \n";
        file << work_lines << std::endl;
        file.close();
    } else {
        throw UnableToWriteToWorkFileException();
    }
}

void update_work_file_postprocessing(const std::string& file_prefix, int n_iter){
    std::string filepath = file_prefix + ".work";

    std::ofstream file(filepath, std::ios::trunc);
    if (file.is_open()) {
        std::string work_lines = "work.program=\"" + carbondate_short_reference() + "\";\n";
        work_lines += "work.operation=\"Post-processing\";\n";
        work_lines += "work.done=100.0; work.passes=" + std::to_string(n_iter) + "; ";
        work_lines += "work.ok=100; \n";
        file << work_lines << std::endl;
        file.close();
    } else {
        throw UnableToWriteToWorkFileException();
    }
}

void check_for_work_file(const std::string& file_prefix) {
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
