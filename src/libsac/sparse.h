#ifndef SPARSEPCM_H
#define SPARSEPCM_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

class SimplePred {
public:
  SimplePred(): lb(0) {}

  double Predict() { return lb; }

  void Update(int32_t val) { lb = val; }

protected:
  int32_t lb;
};

class SparsePCM {
  const double cost_pow = 1;

public:
  SparsePCM(): minval(0), maxval(0), fraction_used(0.0), fraction_cost(0.0) {}

  void Analyse(std::span<const int32_t> buf) {
    // 一次性找到最小值和最大值
    if(buf.empty()) {
      minval = maxval = 0;
      used.clear();
      prefix_sum.clear();
      fraction_used = fraction_cost = 0.0;
      return;
    }

    minval = std::numeric_limits<int32_t>::max();
    maxval = std::numeric_limits<int32_t>::min();
    for(const auto val: buf) {
      if(val < minval) minval = val;
      if(val > maxval) maxval = val;
    }

    // 优化内存布局并统计出现次数
    const size_t range = maxval - minval + 1;
    used.assign(range, false);
    int unique_count = 0;

    for(const auto val: buf) {
      const size_t idx = val - minval;
      if(!used[idx]) {
        used[idx] = true;
        ++unique_count;
      }
    }
    fraction_used = (range > 0) ? (unique_count * 100.0 / range) : 0.0;

    // 预计算前缀和数组
    prefix_sum.resize(range + 1);
    prefix_sum[0] = 0;
    for(size_t i = 0; i < range; ++i) {
      prefix_sum[i + 1] = prefix_sum[i] + (used[i] ? 1 : 0);
    }

    // 优化计算循环
    double sum0 = 0.0, sum1 = 0.0;
    for(const auto val: buf) {
      const int32_t e0 = val;
      const int32_t e1 = map_val(e0);

      sum0 += std::abs(e0); // cost_pow=1 时的优化
      sum1 += std::abs(e1);
    }

    fraction_cost = (sum1 > 0) ? (sum0 / sum1) : 0.0;
  }

  int32_t map_val(int32_t val, int32_t p = 0) const {
    if(val == 0) return 0;

    const int sgn = (0 < val) - (val < 0);
    const int32_t clamped_p = std::clamp(p, minval, maxval);
    const int32_t pidx = clamped_p - minval;

    int32_t start, end;
    if(val > 0) {
      start = pidx + 1;
      end = pidx + val;
    } else {
      start = pidx + val;
      end = pidx - 1;
    }

    start = std::max(start, 0);
    end = std::min(end, static_cast<int32_t>(used.size()) - 1);

    if(start > end) return 0;

    return sgn * (prefix_sum[end + 1] - prefix_sum[start]);
  }

  int32_t minval, maxval;
  double fraction_used, fraction_cost;

protected:
  std::vector<bool> used;      // 使用更紧凑的布尔数组
  std::vector<int> prefix_sum; // 前缀和数组加速区间查询
};

#endif // SPARSEPCM_H