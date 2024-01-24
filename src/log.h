/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_LOG_H
#define CARBONDATE_LOG_H

#include "carbondate_internal.h"


class UnableToWriteToLogFileException : public CarbondateException {
public:
    explicit UnableToWriteToLogFileException(const std::string& file_path) {
        _error_message = "Unable to write to log file at " + file_path;
    }
};

void initialize_log_file();
void update_log_file(const std::string& log_line);
void update_log_file(const std::vector<std::string>& log_lines);

#endif //CARBONDATE_LOG_H
