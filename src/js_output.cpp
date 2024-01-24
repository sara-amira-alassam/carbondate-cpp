/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */
#include <fstream>
#include "js_output.h"


void update_js_output_file(const std::vector<std::string>& output_lines) {
    std::string filepath = project_directory + project_name + ".js";

    std::ofstream file(filepath, std::fstream::app);
    if (!file.is_open()) throw UnableToWriteToOutputFileException(filepath);

    for (const std::string &output_line : output_lines)  file << output_line;

    file.close();
}
