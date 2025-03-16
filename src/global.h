#ifndef GLOBAL_H
#define GLOBAL_H

// #include "windows.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <vector>

#ifndef M_PI
#  define M_PI (3.14159265358979323846)
#endif

// #define NDEBUG

typedef std::vector<double> vec1D;
typedef std::vector<std::vector<double>> vec2D;
typedef std::vector<double*> ptr_vec1D;
typedef std::span<int32_t> span_i32;
typedef std::span<const int32_t> span_ci32;
typedef std::span<const double> span_f64;

#define USE_AVX256
// #define USE_AVX512

#endif
