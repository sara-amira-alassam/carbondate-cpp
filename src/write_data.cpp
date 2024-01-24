/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <fstream>
#include "write_data.h"


std::string log_file_path() {
    return project_directory + project_name + ".log";
}


std::string text_file_path() {
    return project_directory + project_name + ".txt";
}


std::string js_output_file_path() {
    return project_directory + project_name + ".js";
}


void initialize_log_file() {
    std::string filepath = log_file_path();

    update_log_file(carbondate_full_reference());
}


void update_log_file(const std::string& log_line) {
    std::string filepath = log_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToLogFileException(filepath);

    file << log_line << std::endl;
    file.close();
}


void update_log_file(const std::vector<std::string>& log_lines) {
    std::string filepath = log_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToLogFileException(filepath);

    for (const std::string &log_line : log_lines)  file << log_line << std::endl;

    file.close();
}


void update_text_file(const std::vector<std::string>& text_lines) {
    std::string filepath = text_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToTextFileException(filepath);

    for (const std::string &text_line : text_lines)  file << text_line << std::endl;

    file.close();
}


void update_js_output_file(const std::vector<std::string>& output_lines) {
    std::string filepath = js_output_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToOutputFileException(filepath);

    for (const std::string &output_line : output_lines)  file << output_line;

    file.close();
}
