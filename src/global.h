#ifndef GLOBAL_H
#define GLOBAL_H

//#include "windows.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
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

struct SACGlobalCfg {
  static constexpr bool USE_AVX2=true;
  static constexpr int AVX2_MINN=8;
  static constexpr double NLMS_POW_EPS=1.0;
  static constexpr double LMS_ADA_EPS=1E-5;
  static constexpr bool LMS_MIX_INIT=true;// increase stability
  static constexpr bool LMS_MIX_CLAMPW=true;
  static constexpr bool RLS_ALC=true; //adaptive lambda control
};

#endif
