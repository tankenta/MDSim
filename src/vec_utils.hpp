#ifndef VEC_UTILS_H_
#define VEC_UTILS_H_

#include <vector>

const int GRTR = 0;
const int LESR = 1;
const int GEQ = 2;
const int LEQ = 3;
const int EQ = 4;

template <typename T>
std::vector<bool> makeVecMask(std::vector<T> vec, int op, T num);
template <typename T>
std::vector<T> maskedVec(std::vector<T> vec, std::vector<bool> mask);
template <typename T>
std::vector<T> generateRange(T beg, T step, T end);
std::vector<double> linspace(double beg, double end, double n_step);

#endif
