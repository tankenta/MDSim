#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "file_utils.hpp"

void exportPlotData2CSV(
        const std::string& file_name,
        const std::vector<double>& data1, const std::vector<double>& data2,
        const std::string& label1, const std::string& label2) {
    std::ofstream ofs(file_name);
    
    ofs << label1 << "," << label2 << std::endl;
    for (int i = 0; i < data1.size(); i++) {
        ofs << data1[i] << "," << data2[i] << std::endl;
    }
    return;
}
