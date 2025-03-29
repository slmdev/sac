#pragma once

// #include "windows.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <vector>

/*
#ifndef M_PI
#  define M_PI (3.14159265358979323846)
#endif

#define NDEBUG
*/

using vec1D = std::vector<double>;
using vec2D = std::vector<std::vector<double>>;
using ptr_vec1D = std::vector<double*>;
using span_i32 = std::span<int32_t>;
using span_ci32 = std::span<const int32_t>;
using span_f64 = std::span<const double>;

#define USE_AVX256
// #define USE_AVX512
