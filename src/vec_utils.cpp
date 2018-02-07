#include <vector>
#include <cstddef>

#include "vec_utils.hpp"

template <typename T>
std::vector<bool> makeVecMask(std::vector<T> vec, int op, T num) {
    std::vector<bool> mask(vec.size());

    switch (op) {
        case GRTR:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el > num;
            }
            break;
        case LESR:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el < num;
            }
            break;
        case GEQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el >= num;
            }
            break;
        case LEQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el <= num;
            }
            break;
        case EQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el == num;
            }
            break;
    }
    return mask;
}
template std::vector<bool> makeVecMask(std::vector<int> vec, int op, int num);
template std::vector<bool> makeVecMask(std::vector<double> vec, int op, double num);

template <typename T>
std::vector<T> maskedVec(std::vector<T> vec, std::vector<bool> mask) {
    std::vector<T> new_vec;
    const auto head = mask[0];
    for (const auto& el : mask) {
        size_t idx = &el - &head;
        if (el) {
            new_vec.push_back(vec[idx]);
        }
    }
    return new_vec;
}
template std::vector<int> maskedVec(std::vector<int> vec, std::vector<bool> mask);
template std::vector<double> maskedVec(std::vector<double> vec, std::vector<bool> mask);

template <typename T>
std::vector<T> generateRange(T beg, T step, T end) {
    std::vector<T> arr;
    T end_idx = (end - beg) / step;
    for (int i = 0; i <= end_idx; i++) {
        arr.push_back(beg + step*i);
    }
    return arr;
}
template std::vector<int> generateRange(int beg, int step, int end);
template std::vector<double> generateRange(double beg, double step, double end);

std::vector<double> linspace(double beg, double end, double n_step) {
    std::vector<double> arr;
    double step = (end - beg) / (n_step - 1);

    while (beg <= end) {
        arr.push_back(beg);
        beg += step;
    }
    return arr;
}
