#include "work.h"
#include <cstdio>
#include <fstream>

std::string work_file_path() {
    return project_directory + project_name + ".work";
}

void create_work_file() {
    std::ofstream file(work_file_path());
    if (!file.is_open()) {
        throw UnableToCreateWorkFileException(work_file_path());
    }
    file.close();
}

void update_work_file_mcmc(double done, int iterations) {
    std::ofstream file(work_file_path(), std::ios::trunc);
    if (!file.is_open()) throw UnableToWriteToWorkFileException(work_file_path());

    std::string work_lines = "work.program=\"" + carbondate_short_reference() + "\";\n";
    work_lines += "work.operation=\"MCMC\";\n";
    work_lines += "work.done=" + to_string(done * 100., 3) + "; work.passes=" + std::to_string(iterations) + "; ";
    work_lines += "work.ok=100; \n";
    file << work_lines << std::endl;
    file.close();
}

void update_work_file_postprocessing(int n_iter){
    std::ofstream file(work_file_path(), std::ios::trunc);
    if (!file.is_open()) throw UnableToWriteToWorkFileException(work_file_path());

    std::string work_lines = "work.program=\"" + carbondate_short_reference() + "\";\n";
    work_lines += "work.operation=\"Post-processing\";\n";
    work_lines += "work.done=100.0; work.passes=" + std::to_string(n_iter) + "; ";
    work_lines += "work.ok=100; \n";
    file << work_lines << std::endl;
    file.close();
}

void check_for_work_file() {
    std::ifstream file(work_file_path().c_str());
    if (!file.good()) throw WorkFileRemovedException(work_file_path());
}

void remove_work_file() {
    std::remove(work_file_path().c_str());
}
