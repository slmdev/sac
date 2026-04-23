#ifndef RLS_H
#define RLS_H

#include "../global.h"
#include "../common/utils.h"
#include <cmath>

/*
  Model aware adaptive lambda control for RLS
  NIS=err2/(phi+R) with measurement noise R
  for phi->0 reduces to standard norm=err2/E[err2]
*/
class ALC {
  public:
    ALC(double gamma,double beta)
    :gamma(gamma),beta(beta),
    lambda_min(0.99),lambda_max(0.999)
    {
      S0=S1=0.0;
    }
    double Get(double err2,double phi) const
    {
      //unexplained variance is proxy for measurement noise
      double R=std::max(S0-S1,1E-5);
      double nis=err2/(phi+R);//NIS ~ [0,5]
      double m=std::exp(-gamma*nis);
      //double m=1.0-std::tanh(gamma*nis);
      //double m=1.0/(1.0+gamma*nis);
      return lambda_min + (lambda_max-lambda_min)*m;
    }
    void Update(double err2,double phi)
    {
      //update EMAs
      S0 = beta*S0+(1.0-beta)*(err2);
      S1 = beta*S1+(1.0-beta)*(phi);
    }
  private:
    double gamma,beta,lambda_min,lambda_max;
    double S0,S1;
};

// Recursive Least Squares algorithm
class RLS {
  static constexpr double PHI_FLOOR=1E-8;
  public:
    explicit RLS(int n,double gamma,double beta,double nu=1);
    double Predict();
    void Update(double val);
    double Predict(const vec1D &input);
    void UpdateHist(double val);
    int n;
  private:
    double px;
    vec1D x,w;
    vec2D P;
    ALC alc;
};


#endif
