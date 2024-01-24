/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#ifndef CARBONDATE_WRITE_DATA_H
#define CARBONDATE_WRITE_DATA_H

#include "carbondate_internal.h"

class UnableToWriteToLogFileException : public CarbondateException {
public:
    explicit UnableToWriteToLogFileException(const std::string& file_path) {
        _error_message = "Unable to write to log file at " + file_path;
    }
};


class UnableToWriteToTextFileException : public CarbondateException {
public:
    explicit UnableToWriteToTextFileException(const std::string& file_path) {
        _error_message = "Unable to write to text file at " + file_path;
    }
};


class UnableToWriteToOutputFileException : public CarbondateException {
public:
    explicit UnableToWriteToOutputFileException(const std::string& file_path) {
        _error_message = "Unable to write to output file " + file_path;
    }
};


void initialize_log_file();
void update_text_file(const std::vector<std::string>& text_lines);
void update_log_file(const std::string& log_line);
void update_log_file(const std::vector<std::string>& log_lines);
void update_js_output_file(const std::vector<std::string>& output_lines);

#endif //CARBONDATE_WRITE_DATA_H
