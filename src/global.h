#ifndef GLOBAL_H
#define GLOBAL_H

//#include "windows.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <span>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define NDEBUG

using vec1D=std::vector<double>;
using vec2D=std::vector<std::vector<double>>;
using span_i32=std::span<int32_t>;
using span_ci32=std::span<const int32_t>;
using span_cf64=std::span<const double>;

//#define USE_AVX256
//#define UNROLL_AVX256
//#define USE_AVX512

#endif
