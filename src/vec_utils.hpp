#pragma once

#include <cstddef>
#include <vector>

enum class CompOp {
    GRTR, LESS, GEQ, LEQ, EQ
};

template <typename T>
class VecUtils {
public:
    static std::vector<bool> makeVecMask(const std::vector<T>& vec, CompOp c_op, T num);
    static std::vector<T> maskedVec(const std::vector<T>& vec, const std::vector<bool>& mask);
    static std::vector<T> generateRange(T beg, T step, T end);
    static std::vector<double> linspace(double beg, double end, double n_step);
};

template <typename T>
std::vector<bool> VecUtils<T>::makeVecMask(const std::vector<T>& vec, CompOp c_op, T num) {
    std::vector<bool> mask(vec.size());

    switch (c_op) {
        case CompOp::GRTR:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el > num;
            }
            break;
        case CompOp::LESS:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el < num;
            }
            break;
        case CompOp::GEQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el >= num;
            }
            break;
        case CompOp::LEQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el <= num;
            }
            break;
        case CompOp::EQ:
            for (const auto& el : vec) {
                size_t idx = &el - &vec[0];
                mask[idx] = el == num;
            }
            break;
    }
    return mask;
}

template <typename T>
std::vector<T> VecUtils<T>::maskedVec(const std::vector<T>& vec, const std::vector<bool>& mask) {
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

template <typename T>
std::vector<T> VecUtils<T>::generateRange(T beg, T step, T end) {
    std::vector<T> arr;
    T end_idx = (end - beg) / step;
    for (int i = 0; i <= end_idx; i++) {
        arr.push_back(beg + step*i);
    }
    return arr;
}

template <typename T>
std::vector<double> VecUtils<T>::linspace(double beg, double end, double n_step) {
    std::vector<double> arr;
    double step = (end - beg) / (n_step - 1);

    while (beg <= end) {
        arr.push_back(beg);
        beg += step;
    }
    return arr;
}
