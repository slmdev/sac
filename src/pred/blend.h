#ifndef BLEND_H
#define BLEND_H

#include "lms.h"

// blend two expert outputs via sigmoid
template <int log_regression=0>
class Blend2 {
  public:
    Blend2(double beta=0.95,double theta_mu=0.005,double scale=5)
    :scale(scale),mu(theta_mu),
     w(0.5),theta0(0),theta1(1.),
     m0(beta),m1(beta)
    {
    }
    double Predict(double p0,double p1) const
    {
      return w*p0 + (1.0-w)*p1;
    }

    void AdaptTheta(double mu,double delta)
    {
      // pseudo label
      double y = (delta>0)?1.0:0.0;
      double error = y - w;
      double grad0 = error;
      double grad1 = error * delta;

      theta0 += mu*grad0;
      theta1 += mu*grad1;
    }

    void Update(double score0,double score1)
    {
      m0.Update(score0);
      m1.Update(score1);

      double ss0_t=m0.Get();
      double ss1_t=m1.Get();

      // lower score better
      double delta=(ss1_t - ss0_t);

      if constexpr(log_regression) {
        AdaptTheta(mu,delta);
      }

      double z = theta1*delta + theta0;
      z = std::clamp(z,-scale,scale);
      w = 1.0 / (1.0 + std::exp(-z));

   }
  protected:
    double scale,mu;
    double w,theta0,theta1;
    RunSumEMA_NoBC m0,m1;
};

#endif

