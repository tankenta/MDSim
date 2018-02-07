#ifndef FILE_UTILS_H_
#define FILE_UTILS_H_

#include <vector>
#include <string>

void exportPlotData2CSV(
        const std::string& file_name,
        const std::vector<double>& data1, const std::vector<double>& data2,
        const std::string& label1, const std::string& label2);

#endif
