/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#ifndef CARBONDATE_WORK_H
#define CARBONDATE_WORK_H

#include "carbondate_internal.h"


class UnableToCreateWorkFileException : public CarbondateException {
public:
    explicit UnableToCreateWorkFileException(const std::string& file_path) {
        _error_message = "Unable to create writable work file at " + file_path;
    }
};

class UnableToWriteToWorkFileException : public CarbondateException {
public:
    explicit UnableToWriteToWorkFileException(const std::string& file_path) {
        _error_message = "Unable to write to work file at " + file_path;
    }
};

void create_work_file();

void update_work_file_mcmc(double done, int iterations);

void update_work_file_postprocessing(int iterations);

bool work_file_exists();

void remove_work_file();

#endif //CARBONDATE_WORK_H
