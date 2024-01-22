/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_PLAIN_TEXT_H
#define CARBONDATE_PLAIN_TEXT_H

#include "carbondate_internal.h"


class UnableToWriteToTextFileException : public CarbondateException {
public:
    explicit UnableToWriteToTextFileException(const std::string& file_path) {
        _error_message = "Unable to write to text file at " + file_path;
    }
};

void initialize_text_file();
void update_text_file(const std::string& label, const std::vector<double>& values);

#endif //CARBONDATE_PLAIN_TEXT_H
