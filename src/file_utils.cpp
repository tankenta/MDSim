#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "file_utils.hpp"

ProgressBar::ProgressBar(int begin, int end) : begin(begin) {
    size = end - begin + 1;
    beg_time = std::chrono::system_clock::now();
}

void ProgressBar::printProgressBar(int current) {
    int progress = (int)((double)(current - begin)/size * 100);
    auto curr_time = std::chrono::system_clock::now();
    int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(
            curr_time - beg_time).count();
    std::string prog_bar = "";

    for (int i = 0; i < progress/2; i++) {
        prog_bar += "#";
    }
    std::cout << std::setw(3) << progress << "%|"
            << std::left << std::setw(50) << prog_bar
            << std::right << "| ["
            << std::setfill('0') << std::setw(2) << elapsed_sec/60 << ":"
            << std::setw(2) << elapsed_sec%60 << std::setfill(' ') << "]"
            << "\r" << std::flush;
    return;
}

void ProgressBar::finish() {
    auto curr_time = std::chrono::system_clock::now();
    int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(
            curr_time - beg_time).count();
    std::string prog_bar = "";

    for (int i = 0; i < 50; i++) {
        prog_bar += "#";
    }
    std::cout << "100%|" << prog_bar << "| ["
            << std::setfill('0') << std::setw(2) << elapsed_sec/60 << ":"
            << std::setw(2) << elapsed_sec%60 << std::setfill(' ') << "]"
            << std::endl;
    return;
}

// std::chrono::system_clock::time_point printProgressBar(
//         int begin, int end, int current) {
//     int size = end - begin + 1;
//     int progress = (int)((double)(current - begin)/size * 100);
//     static auto beg_time = std::chrono::system_clock::now();
//     auto curr_time = std::chrono::system_clock::now();
//     int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(
//             curr_time - beg_time).count();
//     std::string prog_bar = "";
//
//     for (int i = 0; i < progress/2; i++) {
//         prog_bar += "#";
//     }
//     std::cout << std::setw(3) << progress << "%|"
//             << std::left << std::setw(50) << prog_bar
//             << std::right << "| ["
//             << std::setfill('0') << std::setw(2) << elapsed_sec/60 << ":"
//             << std::setw(2) << elapsed_sec%60 << std::setfill(' ') << "]"
//             << "\r" << std::flush;
//     return beg_time;
// }
//
// void printProgressBar(std::chrono::system_clock::time_point beg_time) {
//     auto curr_time = std::chrono::system_clock::now();
//     int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(
//             curr_time - beg_time).count();
//     std::string prog_bar = "";
//
//     for (int i = 0; i < 50; i++) {
//         prog_bar += "#";
//     }
//     std::cout << "100%|" << prog_bar << "| ["
//             << std::setfill('0') << std::setw(2) << elapsed_sec/60 << ":"
//             << std::setw(2) << elapsed_sec%60 << std::setfill(' ') << "]"
//             << std::endl;
//     return;
// }

void exportPlotData2CSV(
        const std::string& file_name,
        const std::vector<double>& data1, const std::vector<double>& data2,
        const std::string& label1, const std::string& label2) {
    std::ofstream ofs(file_name);
    
    ofs << label1 << "," << label2 << std::endl;
    for (int i = 0; i < data1.size(); i++) {
        ofs << data1[i] << "," << data2[i] << std::endl;
    }

    std::cout << "exported: " << file_name << std::endl;
    return;
}
