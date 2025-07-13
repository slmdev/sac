#ifndef BLEND_H
#define BLEND_H

#include "lms.h"

// blend two expert outputs via sigmoid
class Blend2 {
  public:
    Blend2(double beta=0.95,double theta0=1.0,double theta1=0.0,double scale=5.0)
    :w(0.5),
     th0(theta0),th1(theta1),scale(scale),
     rsum(beta)
    {
    }
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }

    void Update(double score0,double score1)
    {
      rsum.Update(score1-score0);

      double delta=rsum.Get();
      double z = th0*delta + th1;
      z = std::clamp(z,-scale,scale);
      w = 1.0 / (1.0 + std::exp(-z));
   }
  protected:
    double w,th0,th1,scale;
    RunSumEMA rsum;
};

#endif

