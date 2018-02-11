#ifndef FILE_UTILS_H_
#define FILE_UTILS_H_

#include <chrono>
#include <string>
#include <vector>

class ProgressBar {
public:
    ProgressBar(int begin, int end);
    void printProgressBar(int current);
    void finish();

private:
    int begin, size;
    std::chrono::system_clock::time_point beg_time;
};

void exportPlotData2CSV(
        const std::string& file_name,
        const std::vector<double>& data1, const std::vector<double>& data2,
        const std::string& label1, const std::string& label2);

#endif
