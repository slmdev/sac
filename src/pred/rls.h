#ifndef RLS_H
#define RLS_H

#include "../global.h"
#include "../common/utils.h"
#include <cmath>

// adaptive lambda control
template <miscUtils::MapMode tmap_mode>
class ALC
{
  public:
    ALC(double gamma=1.0,double beta=0.95)
    :gamma(gamma),lambda_min(0.99),lambda_max(0.999),
     msum(beta)
    {
    }

    double update(double metric)
    {
      msum.Update(metric);

      // normalize metric by average
      double mnorm = metric/(msum.Get() + 1E-5);
      // map with decay function
      // high mnorm -> low alpha (faster adaption), low mnorm -> high alpha
      double m=miscUtils::decay_map<tmap_mode>(gamma,mnorm);
      return lambda_min + (lambda_max-lambda_min)*m;
    }
  private:
    double gamma,lambda_min,lambda_max;
    RunSumEMA msum;
};

// Recursive Least Squares algorithm
class RLS {
  static constexpr double PHI_FLOOR=1E-8;
  public:
    explicit RLS(int n,double gamma,double nu=1);
    double Predict();
    double Predict(const vec1D &pred);
    void Update(double val);
    void UpdateHist(double val);
    int n;
  private:
    double px,gamma;
    vec1D hist,w;
    vec2D P;
    ALC<miscUtils::MapMode::exp> alc;
};


#endif
