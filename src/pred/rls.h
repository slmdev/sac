#ifndef RLS_H
#define RLS_H

#include "../global.h"
#include "../common/utils.h"
#include <cmath>

// adaptive lambda control
template <miscUtils::MapMode map_mode,int bias_corr=0>
class ALC
{
  public:
    ALC(double gamma=1.0,double beta=0.95)
    :gamma(gamma),beta(beta),power_beta(1.0),
    lambda_min(0.99),lambda_max(0.999),
    eg(0.)
    {
    }

    double update(double metric)
    {
      eg=beta*eg+(1.0-beta)*metric; // EMA

      double eg_hat=eg;
      if constexpr (bias_corr) { // bias correction
        power_beta*=beta;
        eg_hat=eg/(1.0-power_beta);
      };

      // normalize metric by average
      double mnorm = metric/(eg_hat + 1E-5);
      // map with decay function
      // high mnorm -> low alpha (faster adaption), low mnorm -> high alpha
      double m=miscUtils::decay_map<map_mode>(gamma,mnorm);
      return lambda_min + (lambda_max-lambda_min)*m;
    }
  protected:
    double gamma,beta,power_beta;
    double lambda_min,lambda_max;
    double eg;
};

// Recursive Least Squares algorithm
class RLS {
  public:
    explicit RLS(int n,double alpha,double nu=1);
    double Predict();
    double Predict(const vec1D &pred);
    void Update(double val);
    void UpdateHist(double val);
    int n;
  private:
    double px,alpha_mu;
    vec1D hist,w;
    vec2D P;
    ALC <miscUtils::MapMode::exp> alc;
};


#endif
