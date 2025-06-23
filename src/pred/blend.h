#ifndef BLEND_H
#define BLEND_H

#include "lms.h"

// blend two expert outputs via sigmoid
template <int lower_better=1>
class Blend2 {
  public:
    Blend2(double beta=0.95,double alpha=1.0,double scale=5)
    :alpha(alpha),scale(scale),w(0.5),
     m0(beta),m1(beta)
    {}
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }
    void Update(double score0,double score1)
    {
      m0.Update(score0);
      m1.Update(score1);

      double ss0_t=m0.Get();
      double ss1_t=m1.Get();

      double z=0.0;
      if constexpr(lower_better) // lower score better
        z = alpha * (ss1_t - ss0_t);
      else
        z = alpha * (ss0_t - ss1_t);
      z = std::clamp(z,-scale,scale);
      w = 1.0/(1.0+std::exp(-z));
    }
  protected:
    double alpha,scale,w;
    RunSum <> m0,m1;
};

#endif

