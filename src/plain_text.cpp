/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <fstream>
#include "plain_text.h"

std::string text_file_path() {
    return project_directory + project_name + ".txt";
}

void update_text_file(const std::vector<std::string>& text_lines) {
    std::string filepath = text_file_path();

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToTextFileException(filepath);

    for (const std::string &text_line : text_lines)  file << text_line << std::endl;

    file.close();
}
