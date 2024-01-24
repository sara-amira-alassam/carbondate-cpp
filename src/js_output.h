/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_JS_OUTPUT_H
#define CARBONDATE_JS_OUTPUT_H

#include "carbondate_internal.h"

class UnableToWriteToOutputFileException : public CarbondateException {
public:
    explicit UnableToWriteToOutputFileException(const std::string& file_path) {
        _error_message = "Unable to write to output file " + file_path;
    }
};

void update_js_output_file(const std::vector<std::string>& output_lines);

#endif //CARBONDATE_JS_OUTPUT_H
