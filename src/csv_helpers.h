#ifndef CARBONDATE_CSV_HELPERS_H
#define CARBONDATE_CSV_HELPERS_H

std::vector<double> get_csv_data_from_column(
        const std::string& filename,int column_index, char separator);
void write_columns_to_csv(
        const std::string& filename,
        std::vector<std::string> headers,
        std::vector<std::vector<double>> data);
void write_column_to_csv(
        const std::string& filename, const std::string& header, const std::vector<double>& data);
#endif //CARBONDATE_CSV_HELPERS_H
