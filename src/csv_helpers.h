/* carbondate Copyright (C) 2024 Timothy Heaton and Sara Al-Assam
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#ifndef CARBONDATE_CSV_HELPERS_H
#define CARBONDATE_CSV_HELPERS_H

std::vector<double> get_csv_data_from_column(std::fstream* file, int column_index, char separator);
#endif //CARBONDATE_CSV_HELPERS_H
