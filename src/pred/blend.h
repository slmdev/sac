#ifndef BLEND_H
#define BLEND_H

#include "lms.h"

// blend two expert outputs via sigmoid
template <int lower_better=1,int tbias_corr=0>
class Blend2 {
  public:
    Blend2(double beta=0.95,double alpha=1.0,double scale=5)
    :beta(beta),alpha(alpha),scale(scale)
    {
      ss0 = ss1 = 0.0;
      w=0.5;
      power_beta = 1.0;
    }
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }
    void Update(double score0,double score1)
    {
      ss0 = beta * ss0 + (1.0-beta)*score0;
      ss1 = beta * ss1 + (1.0-beta)*score1;

      double ss0_t=ss0;
      double ss1_t=ss1;
      if constexpr(tbias_corr) {
        power_beta *= beta;
        ss0_t = ss0_t / (1.0-power_beta);
        ss1_t = ss1_t / (1.0-power_beta);
      }

      double z=0.0;
      if constexpr(lower_better) // lower score better
        z = alpha * (ss1_t - ss0_t);
      else
        z = alpha * (ss0_t - ss1_t);
      z = std::clamp(z,-scale,scale);
      w = 1.0/(1.0+std::exp(-z));
    }
  protected:
    double w;
    double beta,alpha,scale,power_beta;
    double ss0,ss1;
};

#endif

