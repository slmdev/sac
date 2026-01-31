#ifndef LPC_H
#define LPC_H

#include "../common/utils.h"
#include "../common/math.h"


class OLS {
  public:
    OLS(int n,int kmax=1,double lambda=0.998,double nu=0.001,double beta_sum=0.6,double beta_pow=0.75,double beta_add=2);
    double Predict();
    void Update(double val);
    vec1D x;
  private:
    slmath::Cholesky chol;
    vec1D w,b;
    vec2D mcov;
    int n,kmax,km;
    double lambda,nu,pred;
    double beta_pow,beta_add;
    RunSumGEO esum;
};


#endif // LPC_H
